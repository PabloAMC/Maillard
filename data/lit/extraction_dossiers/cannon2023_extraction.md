# Cannon & Ho 2023 (J. Food Bioact. 21:44-54) — EXTRACTION (file stem `cannon2023`)
### *"Quantitation of key odorants in arginine/cysteine-glucose maillard reactions"*
### Glucose + L-cysteine·HCl (and / or L-arginine), **5 + 5 mmol in 50 mL distilled water (0.1 M each)**, **unbuffered, adjusted to label pH 7.4** with NaOH (Arg-Glu also at its own pH 10.4), 75 mL "reaction vessels" in an oil bath, **1 h at 100 / 130 / 160 °C** (Cys-Glu only at 130 and 160 °C). Volatiles by ether extraction + SAFE, GC-FID with **octanol as internal standard, mg per mol of amino acid, no response factors, no recovery, no replicates stated**. Cys-Glu at 160 °C (final pH 4.4): **MFT 0.43, FFT 0.33, 1-mercaptopropan-2-one 3.97 mg/mol Cys**; all three nd at 130 °C. **Norfuraneol (NF) is not reported in any table.**
DOI: 10.31665/jfb.2023.18338 (printed on p. 44 as 10.31665/JFB.2023.18338)

**Read 2026-11-05 (first reader) from `data/articles/cannon2023.pdf`** (11 pp., publisher PDF, clean text layer).
Prose and Tables 1-4 from `pdftotext -layout`; **Table 3 (the Cys-Glu table) was checked cell by cell against a 130 dpi
render of p. 50** (scratchpad `m4pdf/cannon2023/p-07.png`). Figs. 1-4 are structures and proposed pathways; no figure
carries a quantity. No supplementary material is cited.

**Reservation (M4, case `docs/reference/cases/M4_glucose_mft_protein_amine.md` on branch `m4-case-refresh`, commit
8660c005, sec. 2.9).** The case fixed this paper's role before it was read: a candidate second laboratory for D1 (NF and
MFT in one glucose + Cys pot), else part A evidence (Arg vs Cys). Every thiol value from the Cys-Glu and Arg/Cys-Glu
pots carries **RESERVED: unseen by the engine (M4); no engine run on these pots until M4's design is frozen**
(abbreviated **[RESERVED M4]** below). The role is not met (sec. 9); release is the coordinator's call. The absolute
MFT / FFT levels are also absolute levels from a non-DFA laboratory, so they are additionally marked **[RESERVED H0]**
(H0's tau applies to them; not assessed here). No engine run on any pot of this paper.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[Q]** qualitative statement
by the authors · **[C]** cited by the authors · **[Z]** derived by me, never printed · **[!]** flag. Provenance: **first
reader only, second reader pending.** Nothing here is in a ledger, fit table or benchmark.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Robert J. Cannon*** (McCormick & Company, Hunt Valley, MD; and Dept. Food Science, Rutgers) and **Chi-Tang Ho** (Dept. Food Science, Rutgers University, New Brunswick, NJ) |
| Venue | *Journal of Food Bioactives* **21** (2023) **44-54** (International Society for Nutraceuticals and Functional Foods) |
| DOI | **10.31665/JFB.2023.18338** (printed p. 44), matches the brief |
| Dates | received 20 Feb 2023; revised and accepted 19 Mar 2023 |
| Note | "in commemoration of Professor Michael Granvogl" (author note) |

**Laboratory.** McCormick & Company (Cannon, corresponding) with Rutgers (Ho). Not in `data/keys/laboratories.yml`.
**Proposal, not applied:** a new group `mccormick_cannon_ho`, or, the conservative choice for leave-one-lab-out, merge
into the Ho network (the split files' `ho_network`; `zhang1991b_extraction.md` is the Rutgers Ho glucose + Cys pot).
Owner call. Either way it is **not DFA**: a second laboratory in the sense of the case.

## 1. What this paper answers for the model

| model need | what Cannon 2023 gives | tag |
|---|---|---|
| NF and MFT in one glucose + Cys pot (D1) | **MFT yes (0.43 mg/mol Cys, 160 °C); NF no** (not in any table). No NF / MFT ratio exists | [M] |
| Glucose + Cys MFT at >= 2 pH (D2) | no: one label pH (7.4) | [M] |
| Glucose + Cys MFT vs temperature | 130 °C nd, 160 °C 0.43 mg/mol Cys (1 h); no LOD, so the 130 °C point is censored at an unknown floor | [M] |
| C3 sulfur marker (1-mercaptopropan-2-one, MP) in the same pot | 0.04 (130 °C) / 3.97 (160 °C) mg/mol Cys; MP / MFT ≈ 12 mol/mol at 160 °C [Z] | [M] |
| Amine class at fixed sugar (part A, D7) | Arg-Glu vs Cys-Glu vs Arg/Cys-Glu at 160 °C, label pH 7.4 (sec. 4.4) | [M] |
| Rate constants, time course, pH law | **none**: one time (1 h) per temperature | — |

## 2. Conditions [M] (sec. 2.2, p. 45)

| pot | charge (as printed) | pH as printed | heating |
|---|---|---|---|
| Arg-Glu | glucose 0.90 g (5 mmol) + L-arginine 0.87 g (5 mmol) in distilled water 50 mL | initial 10.4, lowered to **7.4** with 1 M HCl; also run unadjusted at **10.4** (160 °C only) | 1 h at 100 / 130 / 160 °C |
| **Cys-Glu** | glucose 0.90 g (5 mmol) + **L-cysteine hydrochloride 0.79 g (5 mmol)** in 50 mL | initial **1.7**, raised to **7.4** with 2 M NaOH | methods: 1 h at 100 / 130 / 160 °C; **Table 3 and the text report only 130 and 160 °C** [!] |
| Arg/Cys-Glu | glucose 0.90 g + Cys·HCl 0.79 g + L-arginine **0.90 g** ("5 mmol") in 50 mL | initial 5.9, raised to **7.4** with 1 M NaOH | 1 h at 100 / 130 / 160 °C |

| item | as printed | derived [Z] |
|---|---|---|
| Concentrations | 5 mmol each in 50 mL water (+ the NaOH / HCl volume, not stated) | **≈ 0.1 M each** (Cys·HCl 0.79 g = 5.01 mmol anhydrous; glucose 0.90 g = 5.00 mmol; Arg 0.87 g = 4.99 mmol, **0.90 g = 5.17 mmol** in the three-component pot [!]) |
| Buffer | **none** | — |
| pH scale | "initial pH" measured, then adjusted to 7.4: **label pH at room temperature**; after cooling in an ice bath to room temperature "a second pH measurement was made" (**cooled pH**) | — |
| Final (cooled) pH | Cys-Glu 160 °C: **4.4**; Arg/Cys-Glu 160 °C: **4.2**; Arg-Glu 160 °C (label 7.4): **5.2**; Arg-Glu 160 °C (label 10.4): **6.8**. **Not reported for any 100 or 130 °C pot** [!] | Cys-Glu drifts 3.0 pH units over 1 h at 160 °C |
| Vessel / headspace | "75 mL reaction vessels" holding the 50 mL solution; material, closure and atmosphere **not stated** | ≈ 25 mL headspace if one vessel per pot; sealed presumed at 130-160 °C (not stated) [!] |
| Heating / cooling | oil bath, 1 h; then ice bath to room temperature | no come-up time given |
| Replicates | **not stated** | — |

## 3. The assay [M] (secs. 2.3-2.6)

| item | as printed |
|---|---|
| Isolation | whole reaction mixture extracted twice with 25 mL diethyl ether "using an internal standard of octanol at 1 ppm"; Na2SO4; SAFE (< 5 × 10⁻⁴ Torr); N2 to 250 µL |
| GC | Agilent 7890A, FID; apolar RTX-1 50 m and polar CP-Wax 58 (FFAP) 50 m; 40 °C → 2 °C/min → 250 °C |
| Identification | GC-TOF-MS (EI 70 eV, exact mass), against synthesised or purchased standards; retention indices from an ethyl-ester series converted to Kovats by a printed linear equation |
| Quantity reported | **mg per mol** (Table 1 "mg/mol"; Table 3 "mg/mol of cysteine"; Table 4 "mg/mol of arginine") |
| Calibration / response factor | **not stated**: octanol internal standard on FID; no relative response factors, no calibration curves, no SIDA |
| Recovery | **not stated, not measured** |
| LOD / LOQ | **not stated**; "nd: not detected" without a limit |
| Replicate noise | **not stated** (no SD, no n) |

So every value is an internal-standard FID estimate of unknown recovery and response, n unknown. Between-pot ratios of
one compound are on one assay; absolute values are not traceable [!].

## 4. Data

### 4.1 Table 3 (p. 50-51): Cys-Glu, label pH 7.4, 1 h, mg/mol of cysteine [M]

Thiol rows: **[RESERVED M4]**; MFT and FFT also **[RESERVED H0]**.

| compound (as printed) | KI (OV1) | 130 °C | 160 °C |
|---|---:|---:|---:|
| **2-methylfuran-3-thiol (MFT)** | 853 | **nd** | **0.43** |
| **furfuryl mercaptan (FFT)** | 895 | **nd** | **0.33** |
| **1-mercaptopropan-2-one (MP)** | 749 | **0.04** | **3.97** |
| 3-mercaptobutan-2-one | 791 | 0.02 | 9.21 |
| 3-mercaptopentan-2-one | 881 | nd | 0.80 |
| 2-mercaptopentan-3-one | 886 | nd | 1.00 |
| 1-mercaptobutan-3-one / 1-mercaptobutan-2-one | 850 / 852 | nd / nd | 0.36 / 0.76 |
| 2-(1-mercaptoethyl)furan | 932 | nd | 0.47 |
| 5-methyl-2-furanmethanethiol | 994 | nd | 1.21 |
| 3-thiophenethiol | 950 | 0.01 | 5.89 |
| 3-mercapto-2-methylthiophene | 1054 | 0.02 | 2.04 |
| 2-mercaptopropionic acid | 985 | nd | 185.24 |
| 3-mercaptopropionic acid | 998 | 0.01 | 97.38 |
| methyl mercaptan; 1,1-ethanedithiol; 2-mercaptoethanol; 1,2-ethanedithiol | 484-801 | nd | 0.00; 0.14; 0.12; 0.41 |
| 3,5-dimethyl-1,2,4-trithiolane i / ii; 1,2-dithian-4-one | 1112-1124 | nd | 1.31 / 0.60; 2.21 |
| 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone (DMHT) | 1134 | 0.07 | 7.53 |
| 2,5-dimethyl-2,4-dihydroxy-3(2H)-thiophenone (DMDHT) | 1223 | 0.07 | 5.72 |
| thiophene; 2-methyl-; 2,5-dimethyl-; 2-formyl-; 2-acetyl-; 2-formyl-4-methyl-thiophene | 647-1091 | nd | 0.50; 0.16; 0.18; 1.79; 0.33; 4.92 |
| tetrahydrothiophen-3-one; 2-methyltetrahydrothiophen-3-one; thieno[3,2-b]thiophene | 920-1176 | nd | 1.90; 1.02; 1.39 |
| **2-furfural** | 809 | nd | **nd** |
| 2-furfuryl alcohol; 2-acetylfuran; 5-methyl-2-furfural; 5-methyl-2-furanmethanol | 828-947 | nd | 7.14; 1.74; 0.91; 0.66 |
| furaneol (HDMF, "2,5-dimethyl-4-hydroxy-3(2H)-furanone") | 1030 | nd | 2.63 |
| 5-hydroxy-5,6-dihydromaltol | 1116 | 0.14 | 5.20 |
| pyrazine; methyl-; 2,6- + 2,5-dimethyl-; ethyl-; 2,3-dimethylpyrazine | 716-898 | nd | 0.07; 0.94; 0.12; 0.33; 0.16 (total 1.62) |
| 1-hydroxypropan-2-one; acetoin; cyclopentanone; acetaldehyde | 485-763 | nd; nd; nd; 0.08 | 0.29; 0.95; 0.71; 0.06 |
| totals: furans/pyrans; thiophenes; misc. sulfur | — | 0.14; 0.16; 0.07 | 18.28; 33.37; 305.95 |

**Not in Table 3: norfuraneol (4-hydroxy-5-methyl-3(2H)-furanone), 5-HMF (text: "no detected levels"), the
MFT / FFT disulfides.** NF is neither listed as nd nor discussed [!]: the table cannot say whether it was below
detection or not targeted.

### 4.2 Table 4 (p. 52-53): Arg/Cys-Glu, label pH 7.4, 1 h, "mg/mol of arginine" [M] — thiol rows [RESERVED M4]

| compound | 100 °C | 130 °C | 160 °C |
|---|---:|---:|---:|
| **MFT** | **not listed** | **not listed** | **not listed** |
| FFT (furfuryl mercaptan) | nd | nd | 0.30 |
| 1-mercaptopropan-2-one | nd | 0.15 | 3.40 |
| 3-mercaptobutan-2-one | nd | 0.14 | 4.99 |
| 3-thiophenethiol; 3-mercapto-2-methylthiophene | nd | nd | 3.98; 1.74 |
| 2-mercaptopropionic acid; 3-mercaptopropionic acid | nd | nd | 216.46; 42.39 |
| DMHT; DMDHT | nd; 0.00 | nd; 0.00 | 13.57; 36.82 |
| 5-acetyl-2,3-dihydro-1,4-thiazine | nd | 0.07 | nd |
| furaneol (HDMF) | nd | 0.13 | 1.36 |
| 2-furfuryl alcohol; 2-acetylfuran; 5-methyl-2-furfural; 5-methyl-2-furanmethanol | nd | nd | 5.11; 1.85; 1.70; 0.48 |
| 5-hydroxy-5,6-dihydromaltol | nd | 1.15 | 3.52 |
| pyrazine; methylpyrazine | 0.02; nd | 0.02; nd | 0.06; 0.47 |
| other thiophenes (thiophene, 2-methyl-, 2-ethyl-, 2-formyl-, 2-acetyl-, 2-formyl-4-methyl-, tetrahydrothiophen-3-one, 2-methyltetrahydrothiophen-3-one, thieno[3,2-b]thiophene) | nd | nd | 0.67; 0.38; 0.81; 2.51; 2.01; 4.51; 3.22; 1.86; 1.56 |
| totals: furans/pyrans; thiophenes; misc. sulfur | 0.00; 0.00; 0.00 | 1.28; 0.00; 0.36 | 14.03; 73.64; 267.53 |

[!] MFT has no row in Table 4 at all, so "absent" and "not reported" cannot be told apart. [!] "0.00" and "nd" are
both used (DMDHT 0.00 at 100 / 130 °C; total thiophenes 0.00 at 130 °C though DMDHT reads 0.00 there).

### 4.3 Table 1 (p. 47): Arg-Glu, mg/mol, 1 h [M] (part A context; no thiols)

| compound | 100 °C, pH 7.4 | 130 °C, pH 7.4 | 160 °C, pH 7.4 | 160 °C, pH 10.4 |
|---|---:|---:|---:|---:|
| total pyrazines | 0.06 | 0.17 | 7.09 | 207.92 |
| methylpyrazine | nd | 0.01 | 0.42 | 105.40 |
| 2-furfural | nd | 0.49 | 2.96 | 0.63 |
| 2-furfuryl alcohol | nd | 0.10 | 1.57 | 34.65 |
| furaneol (HDMF) | nd | nd | 0.53 | 8.30 |
| 5-hydroxymethyl-2-furfural | nd | 0.04 | 18.58 | nd |
| 1-hydroxypropan-2-one | nd | nd | 0.84 | 16.08 |
| cyclotene | nd | nd | 1.05 | 8.58 |
| 2-aminophenol | nd | nd | 8.11 | 2.70 |
| acetaldehyde | 0.30 | 0.44 | 0.14 | 0.26 |

Remaining Table 1 rows (other pyrazines, quinoxaline, 2-acetylfuran, 5-methylfurfural, 5-methyl-2-furanmethanol,
5-hydroxy-5,6-dihydromaltol, acetoin, cyclohexan-1,2-dione, benzene-1,2-diol) are in the paper; none is an M4 species.
Three hydroxyethylpyrazines (Table 2) are identified by RI and odour only, "< 0.1 mg/mol".

### 4.4 Arithmetic [Z]

1. **Cys-Glu 160 °C, molar** [RESERVED M4; MFT / FFT also RESERVED H0]: MFT 0.43 mg/mol = **3.8 × 10⁻⁶ mol/mol Cys**;
   FFT 2.9 × 10⁻⁶; MP 4.4 × 10⁻⁵; 3-mercaptobutan-2-one 8.8 × 10⁻⁵; HDMF 2.1 × 10⁻⁵ mol/mol.
2. **MP / MFT = 11.7 mol/mol** at 160 °C (within pot, one assay, response factors unknown) [RESERVED M4].
3. **Amine class at 160 °C, label pH 7.4** (per mol of the amino acid named in each table; equimolar charges, Arg 3 %
   high in the three-component pot): adding Arg to Cys-Glu leaves FFT 0.33 → 0.30, MP 3.97 → 3.40, lowers HDMF 2.63 →
   1.36 and 3-mercaptobutan-2-one 9.21 → 4.99, raises DMDHT 5.72 → 36.82 and total thiophenes 33.4 → 73.6; MFT goes
   from 0.43 to not reported. Arg-Glu alone makes 5-HMF 18.6 and 2-furfural 2.96 at 160 °C; with Cys both are nd [M].
   Final pH 4.4 (Cys-Glu) vs 4.2 (Arg/Cys-Glu) vs 5.2 (Arg-Glu): the pots differ in final pH as well as amine.

## 5. What it can and cannot test in the model

**Cannot:**
- D1 (no NF), D2 (one label pH), any rate constant, Ea or pH law (one time per temperature, unbuffered, final pH only at
  160 °C).
- Any absolute comparison without H0's between-laboratory tau, and here the assay has no recovery or response factor.
- 100 °C Cys-Glu: in the methods but absent from the results.

**Can (evidence, after release):**
- A censored temperature contrast for glucose + Cys MFT in one study: nd at 130 °C, 0.43 mg/mol Cys at 160 °C (1 h,
  label pH 7.4, final 4.4), with the floor unknown. Engine status: no engine run on this pot.
- MP in the same pot, rising ~100× from 130 to 160 °C, as a C3-sulfur marker beside MFT (rival R3's intermediate).
- Part A: Arg vs Cys vs Arg + Cys at fixed glucose and label pH (sec. 4.4 item 3), one laboratory, final pH confounded.

## 6. Conflicts and agreements with existing records

- **Hofmann & Schieberle 1998 Table 2** (`hofmann1998_reconciliation.md` sec. 2.2; DFA, SIDA; Cys 33 + glucose 100 mM,
  0.5 M phosphate, 145 °C, 20 min): glucose MFT 0.3 / 1.9 / 0.4 µg per 3.3 mmol Cys at pH 3 / 5 / 7, i.e. 0.09 / 0.58 /
  0.12 mg/mol Cys [Z]. Cannon's 0.43 mg/mol Cys (160 °C, 1 h, label 7.4 → 4.4, 1 : 1 Cys : glucose, no buffer) is the
  same order **[RESERVED M4] [RESERVED H0]**: an across-laboratory absolute comparison, recorded, not scored. No conflict
  claimed.
- **`zhang1991b_extraction.md`** (Rutgers Ho, glucose + Cys 50 / 50 mM, label 7.5, 180 °C, SDE): no MFT or FFT detected.
  Cannon (same network, SAFE, 160 °C) detects both. The authors note their total pyrazines (1.6 mg/mol) are far below
  Zhang & Ho's 117.6 mg/mol and blame SDE vs SAFE [Q]. Not a conflict in kind: different T, extraction and no LODs.
- **`deroos2005_extraction.md`** (Givaudan, 1 : 1 Cys + glucose, lactate pH 5, ~100 °C): glucose MFT < 0.05 mg/mol Cys.
  Consistent with Cannon's nd at 130 °C (no LOD there). Cannon cites de Roos for low hexose thiols [C].
- **[!] Internal:** Cys-Glu at 100 °C promised in the methods, absent from Table 3 and the text; Arg 0.90 g labelled
  "5 mmol" (5.17 mmol); "nd" vs "0.00"; Table 4's unit is "mg/mol of arginine" although the thiols come from cysteine
  (numerically ~equal at equimolar charge).
- Compound keys: MFT (`2_methyl_3_furanthiol`), FFT (`2_furfurylthiol`) resolve in `data/keys/compounds.yml`; 1-mercapto-
  2-propanone not checked; nothing was keyed (no data rows).

## 7. Provenance tiers (first reader only)

All table cells: `direct_measurement`, quant class **semi-quantitative (internal standard on FID, no response factor, no
recovery, n not stated)**, confidence low. Every [Z] value: `derived`. **Proposed ledger entries: none** until the
reservation is released; if released, at most one sulfur-ledger row for the Cys-Glu 130 / 160 °C thiols with flags
`no_recovery, no_response_factor, no_lod, n_not_stated, unbuffered, final_ph_only_at_160c`.

## 8. Held dossiers that already cite this paper

None. `docs/guides/CITATION_CONTEXT_2026-10.md` lists it as MAYBE ("an abstract that says quantified"), not read until now.

## 9. M4 role (case 8660c005)

**Role fixed before reading** (case sec. 2.9, candidate table): "decisive for D1 (second laboratory) only if NF and
MFT are both quantified in one Cys + glucose pot; else part A evidence (Arg vs Cys)". The measurement that would count:
"NF and MFT in one pot, one assay each with stated recovery".

| element of the criterion | what the paper has | verdict |
|---|---|---|
| one Cys + glucose pot | yes: 0.1 M + 0.1 M, water, label pH 7.4, 1 h at 130 and 160 °C | met |
| MFT quantified in it | yes at 160 °C: **0.43 mg/mol Cys**; nd at 130 °C (no LOD) | met (160 °C only) |
| NF quantified in the same pot | **no: norfuraneol does not appear in Table 3 or the text** | **not met** |
| one assay each with stated recovery | **no: octanol internal standard on FID; no recovery, no response factor, no LOD, n not stated** | **not met** |

**Verdict: NOT MET** (the criterion requires NF and MFT both; NF is absent, and no recovery is stated). D1 keeps one
laboratory (DFA, Hofmann 1998 Table 5). The fallback role applies: **part A evidence**, Arg vs Cys vs Arg + Cys at fixed
glucose (sec. 4.4 item 3), one laboratory, final pH confounded.

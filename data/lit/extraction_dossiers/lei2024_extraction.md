# Lei, Wang, Zhao, Dou, Zhang, Wang, Gao, Wang, Xu & Dong 2024 (Food Res. Int. 188:114454) — EXTRACTION (file stem `lei2024`)
**DOI:** 10.1016/j.foodres.2024.114454 (printed on p. 1; the PDF metadata carry the same DOI and title)
### *"Analysis of volatiles from the thermal decomposition of Amadori rearrangement products in the cysteine-glucose Maillard reaction and density functional theory study"*
### A purified cysteine-glucose "ARP" (0.01 mmol, solid, no water) heated in a micro-chamber thermal extractor at 120 - 180 C for 120 min. Volatiles were purge-trapped on Tenax, then analysed by TD-GC-MS, semi-quantified against one internal standard (RF 1). **MFT appears at 150 and 160 C only (4.00 ± 1.11 and 2.04 ± 0.03 ug/g), identified by RI and MS only, with an RI (979) 110 units from the literature value it cites (869). No norfuraneol. No recovery. No glucose + Cys comparison pot.** The second half is a DFT study (excluded by the no-computed-barrier rule), and the SI holds only DFT rate constants.

**Read 2026-11-06 (first reader)** from `lei2024.pdf` (12 pp., publisher PDF, clean text layer) and
`lei2024_supplementary.docx` (one table). Text from `pdftotext -layout` and `textutil`. Table 1 was read from the text
layer; its column alignment was checked against the printed column order. No figure was digitised: Figs. 2 - 3 are a
chromatogram, a heatmap and PCA plots, and Figs. 4 - 6 are DFT schemes.

**Role, fixed before reading** (case M4, `docs/reference/cases/M4_glucose_mft_protein_amine.md` sec. 2.9, sweep 2's
first candidate; confirmed in commit a29535db before this file was opened):
- **decisive for D1 in its fed form only if** NF and MFT are both quantified from the fed Cys-glucose ARP with stated
  recovery, and the laboratory is not one whose pots are M2's TEST (BTBU, CAAS);
- **otherwise evidence** for D5 and the arm-3-type question (MFT per mol of fed Amadori; R4' vs the claim);
- any computed (DFT) part is excluded.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown · **[Z]**
derived by me · **[!]** flag · **[SI]** from the supplementary material. Provenance: **first reader only, second reader
pending.** Nothing here is in a ledger, fit table or benchmark. No engine run.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Liming Lei, Shang Wang (equal contribution), Ziqian Zhao, **Shaohua Dou** (corresponding), Sufang Zhang, Yi Wang, Pengxun Gao, Binchen Wang, Xianbing Xu, **Liang Dong** (corresponding) |
| Affiliations | School of Food Science and Technology, Dalian Polytechnic University (National Engineering Research Centre of Seafood), Dalian; College of Life and Health, Dalian University; Dalian Huamei School |
| Venue | *Food Research International* **188** (2024) 114454 |
| DOI | 10.1016/j.foodres.2024.114454 |
| Dates | received 26 Jan 2024; revised 28 Apr 2024; accepted 1 May 2024; online 3 May 2024 |
| Supplementary | `lei2024_supplementary.docx`: "Table SM", DFT rate constants only |

**Laboratory.** Dalian Polytechnic University, Dong / Xu group (the authors' earlier DFT papers on lipid oxidation). It
is **not** BTBU or CAAS (not one of D3's laboratories), and it is not in `data/keys/laboratories.yml`. Proposed id
`dalian_dpu_dong`.

## 1. What this paper answers for the model

| model need | what Lei 2024 gives | tag |
|---|---|---|
| NF and MFT from a fed Cys-glucose ARP (D1, fed form) | **MFT only** (150 / 160 C), **no norfuraneol** in the 60 volatiles of Table 1 | [M] |
| MFT per mol of fed Amadori consumed against a matched glucose + Cys pot (arm 3 / R4') | **no**: ARP consumption not measured, no sugar + Cys pot | -- |
| the ARP's identity and purity | NMR shifts listed, with ambiguous integrals ("0H") and a 17.37 ppm carbon that no Fru-Cys structure carries [!] | [M] / [!] |
| temperature dependence of the volatiles | Table 1, seven temperatures (with a 100 C column), 120 min each | [M] |
| rate constants | **DFT only** (Fig. 6, Table SM): excluded | -- |

## 2. Conditions [M]

| item | as printed |
|---|---|
| ARP synthesis | glucose 0.2 mol + NaHSO3 2 g in ethanol 60 mL / glycerol 30 mL; Cys 0.06 mol + acetic acid 8 mL added at 105 C, 5 h at 105 C (Glinsky 1996 method [C]); Dowex 50WX4-H+ (water eluate until TTC-negative, then 1 M ammonia); rotary evaporation, freeze-drying; semi-preparative C18 HPLC (5 - 100 % methanol), UV 271 nm |
| heating pot | **0.01 mmol ARP** (solid) in 20 mL glassware with 10 uL o-dichlorobenzene (50 mg/L) as internal standard, in an M-CTE250 micro-chamber thermal extractor; three chambers in parallel |
| temperature / time | 120, 130, 140, 150, 160, 170, 180 C for 120 min each (Table 1 also prints a 100 C column) |
| medium / water / pH | **none stated**: a dry solid; no water, no buffer, no pH |
| atmosphere | air: chamber closed during heating, then purged with air 100 mL/min for 30 min onto Tenax TA |
| replicates | "All experiments were conducted in triplicate" [V]; mean ± SD (Table 1 footnote says "standard error") [!] |

## 3. Assay [M]

- TD (280 C, 5 min) -> cold trap (10 C) -> GC-MS (Agilent 6890 / 5975C, HP-5ms, splitless, EI 70 eV, m/z 40 - 400).
- **Quantification:** "relative content ... according to the proportion of peak area and the content of internal
  standard" [V]. That is one internal standard (o-dichlorobenzene) with response factor 1 for every compound: a
  **semi-quantitative** basis. Units ug/g (of ARP, by context [Z]).
- **No recovery, LOD or calibration** for any volatile.
- Identification: "RI, MS, STD" for most compounds. **MFT (no. 46): "RI, MS" only**, with no standard.

## 4. Data

### 4.1 Table 1 (pp. 4 - 6), selected rows, ug/g, mean ± SD (n = 3) [M]

| no. | compound | RI / LRI (cited) | 100 C | 120 C | 130 C | 140 C | 150 C | 160 C | 170 C | 180 C |
|---|---|---|---|---|---|---|---|---|---|---|
| 46 | **2-methyl-3-furanthiol** | **979 / 869** (Majcher & Jelen 2007) | ND | ND | ND | ND | **4.00 ± 1.11** | **2.04 ± 0.03** | ND | ND |
| 26 | 2-methylthiophene | 638 / 778 | ND | 2.85 ± 0.29 | ND | 12.44 ± 5.80 | 27.84 ± 8.14 | ND | 9.90 ± 1.20 | ND |
| 13 | 2-acetylthiazole | 1018 / 1011 | 0.65 ± 0.15 | 3.87 ± 0.55 | 3.09 ± 0.67 | 25.98 ± 20.87 | ND | 8.36 ± 0.65 | 16.04 ± 2.57 | 17.19 ± 0.35 |
| 43 | 2-ethylfuran | 617 / 691 | 4.24 ± 1.23 | 6.80 ± 1.50 | 4.48 ± 1.2 | 13.65 ± 2.36 | 32.46 ± 6.14 | 9.72 ± 0.76 | 9.26 ± 0.54 | 11.99 ± 1.84 |
| -- | norfuraneol, furfural, FFT | -- | not in Table 1 at any temperature | | | | | | | |

Text: total sulfur volatiles 63.05 ug/g at 120 C to 555.95 at 180 C [M]. "Sulphur-containing heterocyclic volatiles
could not be detected at temperatures below 120 C" in earlier experiments (data not given) [Q]. That is inconsistent
with the 100 C column's 2-acetylthiazole [!].

### 4.2 Arithmetic [Z]

- MFT 4.00 ug/g of ARP at 150 C, on Fru-Cys (C9H17NO7S, 283.3 g/mol), is ~9.9e-6 mol MFT per mol ARP charged.
  Hofmann 1998's glucose + Cys pot (1.9 ug from 3.3 mmol Cys, 145 C, 20 min, aqueous) gives ~5.0e-6 mol per mol Cys.
  Same order of magnitude, but different matrix (dry vs 0.5 M phosphate), time, atmosphere and quantification
  (RF 1 vs SIDA), so the two are not comparable as an arm-3 ratio.
- The MFT RI here (979) is 110 units above the cited LRI (869) and above Ames 2001's 873 on HP-5. Other rows show
  similar or larger RI-LRI gaps (2-methylthiophene 638 / 778; pyrrole 634 / 1508). So the RI column is not a reliable
  identity check in this paper, and MFT's identity rests on a library MS match [!].

### 4.3 Supplementary material [SI]

"Table SM: Variation of reaction rate constants with temperature" [SI]. DFT rate constants (s-1) for the 2-methylthiophene
pathway steps at 273 - 453 K. **Excluded** by the project's no-computed-barrier rule (`assert_no_dft_*` guards; AGENTS.md).
It has no measured value.

## 5. What it can and cannot test in the model

**Can (weak evidence only):**
- A purified Cys-glucose Amadori-type compound gives sulfur heterocycles from 120 C, and a compound library-matched to
  MFT at 150 - 160 C, dry, under air.
- Qualitatively, a fed hexose-cysteine pool can feed sulfur volatiles at process temperature.

**Cannot:**
- D1 (no NF);
- arm 3's question (no consumption of the ARP, no matched sugar pot);
- any rate or yield on a quantitative basis (RF 1, no recovery);
- any aqueous condition, pH or buffer.

Also flagged:
- The ARP's identity is not secured. Cysteine + glucose is known to give 2-(polyhydroxyalkyl)thiazolidine-4-carboxylic
  acids as well as the Amadori compound; the 17.37 ppm carbon and the "0H" integrals are unexplained [!].
- The DFT mechanism (Fig. 4: 3-DG under acid, "1-deoxyglucuronide" [sic, 1-DG] under neutral, then glyceraldehyde +
  acetaldehyde + H2S -> 2-methylthiophene) is a computed proposal and is not used.

## 6. Conflicts and agreements with existing records

1. **Agrees with de Roos 2005 Table II and Du 2024b** (`du2024b_extraction.md`): a Cys-hexose Amadori-type compound is a
   real, isolable pool.
2. **Bears on case M4's D5 / R4'** only qualitatively. A fed pool gives a compound library-matched to MFT, which is
   what R1 / R5 expect. R4' (the pool as a Cys sink) also releases Cys and sugar on decay, so a fed-pool MFT signal
   without a matched sugar pot and a consumed-ARP denominator does not separate them.
3. [!] RI / LRI mismatches throughout Table 1; the SD vs SE wording; the 100 C column against the text.

## 7. Provenance tiers (first reader only; second reader pending)

- Table 1 values: `semi_quantitative_single_IS` (RF 1, no recovery), confidence **low** for MFT (RI / MS only, RI
  mismatch), low-medium for compounds with "STD".
- ARP identity: `asserted_by_nmr`, confidence low.
- SI: computed, excluded.
- **Proposed ledger entries: none.** No row can be simulated: dry solid, air purge, no water basis.

## M4 role (case M4 sec. 2.9; roles confirmed before reading, commit a29535db)

| element of the criterion | what the paper has | verdict |
|---|---|---|
| a fed Cys-glucose ARP pot | yes (purified, identity not secured) | met with a flag |
| MFT quantified from it | at 150 and 160 C only; RI / MS only, RF 1 | met weakly |
| **NF quantified in the same pot** | **no** (norfuraneol not in Table 1) | **not met** |
| stated recovery | **no** | **not met** |
| laboratory not one of M2's TEST labs (BTBU, CAAS) | yes (Dalian Polytechnic) | met |

**Verdict: decisive role NOT MET.** The fallback role applies: evidence for D5 (a fed pool gives sulfur volatiles and a
compound library-matched to MFT at 150 - 160 C, dry). It does not answer the arm-3 question: no consumed-ARP
denominator, no matched glucose + Cys pot. It does not separate R4' from R1 / R5.

It changes nothing in the case's status ("ready to freeze as development grade, decider W6") or in the plan to build
only the W3 structure plus the reduced-model power check.

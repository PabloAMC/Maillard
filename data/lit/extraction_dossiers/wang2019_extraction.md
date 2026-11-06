# Wang, Hu, McClements, Nie, Shen, Li, Huang, Zhong, Chen, Zeng & Xie 2019 (Food Res. Int. 123:403-413) — EXTRACTION (file stem `wang2019`)
**DOI:** 10.1016/j.foodres.2019.05.001 (printed on p. 403; the PDF metadata carry the same DOI and title)
### *"pH and lipid unsaturation impact the formation of acrylamide and 5-hydroxymethylfurfural in model system at frying temperature"*
### Asn 50 mM + glucose 25 mM, 2 mL, in **unbuffered water (initial pH 4.6) or 0.4 M sodium acetate (pH 6.0)**, with or without a layer of free fatty acid or triglyceride on top, in sealed pressure tubes in an oven at **180 C for 30 min**, n = 3.
### Measured: acrylamide (SIDA, 13C3), HMF (external standard), fructose-asparagine and the alpha-dicarbonyls 3-DG / GO / MGO / DA (OPD, external quinoxaline standards); cooled pH after heating (SI Table S7).
### **The lipid-free control (AG-control) gives HMF at two pH with lipid held fixed: ~12 vs ~9 nmol/mL (Fig. 2C, figure only). But pH and buffer change together (water vs 0.4 M acetate), and the cooled pH converges (5.7 vs 5.9).**

**Read 2026-11-06 (first reader)** from `wang2019.pdf` (11 pp., publisher PDF, text layer clean except figures) and
`wang2019_supplementary.docx` (Tables S1 - S8, Figs. S1 - S3). Text from `pdftotext -layout` and `textutil`. Tables 2
and 3 were read from the text layer and checked on a 100 dpi render of p. 410. **Figs. 1C and 2C (the AG-control bars)
were read by eye** from 100 - 110 dpi renders of pp. 406 - 407: read error about ±1 nmol/mL for HMF and ±10 nmol/mL for
acrylamide. Fig. 3 (Fru-Asn) was not digitised.

**Role, fixed before reading** (M0a's sweep candidate, carried into case M4 part A,
`docs/reference/cases/M4_glucose_mft_protein_amine.md` sec. 1; confirmed in commit a29535db before this file was opened):
- **decisive for part A (HMF) only if** an Asn + hexose pot carries HMF at >= 2 pH with the lipid held fixed;
- otherwise evidence.

Part A is an engine-structure requirement only (FM7), so a met role feeds a later part-A case, not M4.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[Z]** derived by me · **[!]** flag · **[SI]** from the
supplementary material. Provenance: **first reader only, second reader pending.** Nothing here is in a ledger, fit
table or benchmark. No engine run.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Yuting Wang, Huiyu Hu, David Julian McClements, Shaoping Nie, **Mingyue Shen** (corresponding), Chang Li, Yousheng Huang, Yadong Zhong, Jie Chen, Maomao Zeng, **Mingyong Xie** (corresponding) |
| Affiliations | State Key Laboratory of Food Science and Technology, Nanchang University; Dept. of Food Science, University of Massachusetts Amherst; State Key Laboratory of Food Science and Technology, Jiangnan University; Jiangxi Institute of Analysis and Testing |
| Venue | *Food Research International* **123** (2019) 403-413 |
| DOI | 10.1016/j.foodres.2019.05.001 |
| Dates | received 28 Nov 2018; revised 12 Apr 2019; accepted 1 May 2019; online 2 May 2019 |
| Supplementary | `wang2019_supplementary.docx`: Tables S1 - S6 (MS parameters), S7 (cooled pH), S8 (two-factor ANOVA); Figs. S1 - S3 (chromatograms) |

**Laboratory.** Nanchang University, Xie / Shen group (State Key Laboratory, Nanchang); the experiments were run
there. Not in `data/keys/laboratories.yml`. Proposed id `nanchang_xie`. Independent of the Schibilsky / Kroh / Kanzler
laboratories behind M0a's spent pots.

## 1. What this paper answers for the model

| model need | what Wang 2019 gives | tag |
|---|---|---|
| HMF from Asn + glucose at >= 2 pH, lipid fixed (part A's role) | AG-control: unbuffered (initial pH 4.6) ~12 nmol/mL vs 0.4 M acetate (pH 6.0) ~9 nmol/mL, 180 C, 30 min (Fig. 2C) | [F] |
| the confound | the medium changes with pH (water vs 0.4 M acetate); the cooled pH after heating is 5.7 ± 0.1 vs 5.9 ± 0.04 (Table S7) | [M] [SI] |
| alpha-dicarbonyls in the same control pots | 3-DG 65.6 vs 24.6; GO 22.6 vs 23.79; **MGO 65.7 vs 322**; DA 54.88 vs 67.1 nmol/mL (Table 2) | [M] |
| acrylamide in the same control pots | ~228 vs ~278 nmol/mL (Fig. 1C) | [F] |
| lipid effects (M5) | Figs. 1 - 3, Tables 2 - 3 | [M] / [F] |
| thiol, MFT, Cys | **none** | -- |

## 2. Conditions [M]

| item | as printed |
|---|---|
| charge | Asn : Glc 2 : 1, **50 mM : 25 mM**, 2 mL in 15 mL thick-wall pressure tubes (PTFE screw cap); lipids layered on top (50 / 100 / 150 umol FFA, or 50 umol TG); the AG-control has none |
| medium / pH | **non-buffered: distilled water, pH 4.6** (attributed to asparagine) or **buffered: 0.4 mol/L sodium acetate-acetic acid, pH 6.0**; label pH, scale (room temperature) not stated for the initial values |
| cooled pH after heating | Table S7 [SI]: AG-control **5.7 ± 0.1 (non-buffered), 5.9 ± 0.04 (buffered)**; with lipids 4.6 - 5.2 (non-buffered), 5.8 - 5.9 (buffered); measured at room temperature |
| heating | electric blast oven, **180 ± 2 C, 30 min**; ice-bath cooling; vessel heat-up not stated |
| headspace / O2 | air (sealed tubes, ~13 mL headspace over 2 mL [Z]); not discussed |
| replicates | three independent replications; Table 2 footnote "mean ± standard deviation (n = 3)", Methods "mean ± standard error" [!] |

## 3. Assays [M]

- **Acrylamide:** UPLC-QqQ-MS/MS, 13C3-acrylamide internal standard (50 - 1000 ng/mL). A SIDA basis.
- **HMF:** same run, **external standard** (5 - 100 ng/mL, r2 0.998). No internal standard, no recovery stated.
- **Fru-Asn:** UPLC-MS/MS, external standard (10 - 800 ng/mL).
- **alpha-Dicarbonyls:** post-incubation OPD (50 uL sample + phosphate pH 7 + OPD, 25 C, 3 h, dark), HPLC-QqQ-MS/MS;
  external quinoxaline calibrations (3-DG 2 - 30, GO 0.5 - 12, MGO 5 - 120, DA 0.5 - 12 nmol/mL).
- MS parameters: Tables S1 - S6 [SI]; chromatograms Figs. S1 - S3 [SI].
- **No LOD / LOQ** is printed for any analyte.
- [!] The buffered control's MGO (322 nmol/mL) is above the printed MGO calibration top (120 nmol/mL); dilution is
  not stated.

## 4. Data

### 4.1 The lipid-free control, both media (180 C, 30 min)

| analyte | non-buffered (initial pH 4.6; cooled 5.7) | buffered 0.4 M acetate (pH 6.0; cooled 5.9) | ratio non-buf / buf | source |
|---|---|---|---|---|
| HMF | ~12 ± 1 nmol/mL | ~9.2 ± 0.3 | ~1.3 (+0.12 dex) | Fig. 2A / B / C [F]; marked significant (*) |
| acrylamide | ~228 ± 10 | ~278 ± 3 | ~0.82 (-0.09 dex) | Fig. 1C [F]; significant (*) |
| 3-DG | 65.6 ± 4.3 | 24.6 ± 1.8 | 2.67 (+0.43 dex) | Table 2 [M] |
| GO | 22.6 ± 2.0 | 23.79 ± 0.89 | 0.95 | Table 2 [M] |
| MGO | 65.7 ± 5.5 | 322 ± 15 | 0.20 (-0.69 dex) | Table 2 [M] |
| DA | 54.88 ± 0.81 | 67.1 ± 1.0 | 0.82 | Table 2 [M] |
| Fru-Asn | higher in non-buffered (Fig. 3) | lower | -- | [Q], not digitised |

The SI's Table S8 two-factor ANOVA gives p < 0.001 for pH, lipid type and their interaction on every analyte [SI].

### 4.2 Arithmetic [Z]

- **Per unit of initial pH** (4.6 -> 6.0): HMF -0.09 dex/pH, 3-DG -0.30 dex/pH, MGO +0.49 dex/pH.
- **Per unit of cooled pH** (5.7 -> 5.9): the same contrasts are ~7x steeper per unit. **The pH axis is not
  identified**, because the two arms differ in buffer (0.4 M acetate vs none) as much as in pH, and their end pH nearly
  converge.
- The MGO rise with pH (x4.9) has the sign of Yu 2017's MGO vs pH (0.50 - 0.76 dex/pH, 120 C). The 3-DG fall with pH
  has the sign of acid-favoured 3-DG.

## 5. What it can and cannot test in the model

**Can (evidence; a later part-A case, not M4):**
- the acrylamide lane's HMF, 3-DG, MGO and acrylamide at two media at 180 C in an Asn + glucose pot. These are
  quantities the engine can represent (acrylamide lane + trunk dicarbonyls) once part A's one pH module exists;
- a buffer x pH contrast for part A's pH module, if the buffer term (acetate catalysis, M5) is declared.

**Cannot:**
- a pH law on its own (pH and buffer are confounded; the end pH converges);
- anything on thiols, cysteine or MFT (no Cys);
- kinetics (one time point).

## 6. Conflicts and agreements with existing records

1. **Agrees in sign with Yu 2017** (`yu2017_extraction.md`): MGO rises with pH in a glucose + amino-acid pot.
2. **Agrees with the acrylamide lane's own pH factor in sign** (De Vleeschouwer 2006, B15): acrylamide higher at the
   higher pH. Here it is confounded with 0.4 M acetate.
3. **M0a / CF-133 (path P1):** in this pot the engine's acrylamide lane applies no trunk pH term to the Amadori steps,
   and Asn's Amadori chain (Fru-Asn, measured here) is exactly the chain P1 / P1b concern.
4. [!] MGO above the calibration top; SD vs SE wording.

## 7. Provenance tiers (first reader only; second reader pending)

- Table 2 / 3 values: `direct_measurement`, external quinoxaline calibration. Confidence medium (the buffered MGO is
  above the calibration range: low).
- Fig. 1C / 2C control bars: `figure_read_by_eye`, ±1 (HMF), ±10 (acrylamide) nmol/mL. Confidence medium-low.
- Table S7 pH: `direct_measurement`, cooled, room temperature. Confidence medium.
- **Proposed ledger entry (not written):** `wang2019_asn_glc_control_two_media_180c`. Asn 50 + Glc 25 mM, water vs
  0.4 M acetate, 180 C, 30 min. Species HMF, acrylamide, 3-DG, GO, MGO, DA. Flags: ph_and_buffer_confounded,
  ph_cooled_reported, hmf_external_standard, mgo_above_calibration. Laboratory `nanchang_xie`.

## M4 role (part A, inherited from M0a; roles confirmed before reading, commit a29535db)

| element of the criterion | what the paper has | verdict |
|---|---|---|
| an Asn + hexose pot | yes: Asn 50 + glucose 25 mM | met |
| HMF at >= 2 pH | yes: initial pH 4.6 and 6.0 (cooled 5.7 and 5.9) | met in form |
| lipid held fixed | yes: the AG-control has no lipid in either medium | met |
| one design, pH the only factor | **no**: water vs 0.4 M acetate | **confounded** |

**Verdict: the decisive role is MET IN FORM but CONFOUNDED.** It can be part A's HMF pH datum only with:
- a declared acetate-buffer term (M5's buffer factor);
- the pH scale stated (initial label vs cooled), since the two arms end 0.2 units apart.

For M4 it changes nothing: part A is structure only (FM7). It is recorded for the later part-A case, as evidence with
this confound.

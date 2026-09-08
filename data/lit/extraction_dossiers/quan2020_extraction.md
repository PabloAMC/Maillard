# Quan et al. 2020 — EXTRACTION (lysine, asparagine, tryptophan and their mixture + glucose in 0.1 M phosphate pH 7.0, 100 and 130 C, 3-21 min; glyoxal, methylglyoxal, CML, CEL, acrylamide, harmane, norharmane, melanoidins; eleven fitted constants without units)
### The only sugar + free-amine pot on disk with glyoxal AND methylglyoxal measured at 100 and 130 C in water — and a rate-constant table the repository has already refused; this dossier records what survives of it.

**Source on disk:** `data/articles/quan2020.pdf` (owner's download, file dated 2026-08-28, read 2026-09-08;
Elsevier journal pre-proof, 40 pages, Food Chemistry, PII S0308-8146(20)31249-8, DOI
10.1016/j.foodchem.2020.127387; version of record Food Chem. 332 (2020) 127387 — NOT on disk). Read from
the `pdftotext` text layer in the scratchpad; Table 1 (pp. 32-33) re-extracted with pypdf layout mode to fix
the column positions — the dash cells then align, and the assignment of the continuation columns below rests
on that alignment plus the text. Figures 1-5 (pp. 34-38) are images: all concentration-time data and the
reaction scheme are **FIGURE-ONLY**. Supplementary Figures 1 (rate equations) and 2, and Supplementary
Tables 1-2 (validation), are not on disk. Repo status before this dossier: **REFUSED** for parameters in
`k1_kinetic_parameters.md` §2d, `k3_final_parameter_inventory.md` C.5 and row "Quan 2020 Table 1, all 88
constants", and `docs/reference/FIT_HOLDOUT_DECLARATION.md` line 68 ("neither — no units, no orders");
`parameters_acrylamide.py` line 720 lists "Quan 2020, all 88 rate constants" as refused. This dossier does
not overturn that verdict; it adds the levels B18 section 6 and Programme 7 asked for.

## 0. Identity

| field | value |
|---|---|
| Title | "Simultaneous generation of acrylamide, β-carboline heterocyclic amines and advanced glycation end products in an aqueous Maillard reaction model system" |
| Authors | Wei Quan, Yong Li, Ye Jiao, Chaoyi Xue, Guoping Liu, Zhaojun Wang, Zhiyong He, Fan Qin, Maomao Zeng (corresponding), Jie Chen (corresponding); Jiangnan University, Wuxi |
| Venue | Food Chemistry, accepted 17 Jun 2020 (received 28 Mar, revised 14 Jun 2020); journal pre-proof |
| Systems | Lys/Glu, Asn/Glu, Trp/Glu, Mix/Glu (Lys + Asn + Trp + Glu) |
| Naming | "MRHPs" = Maillard reaction harmful products; "HAs" = β-carboline heterocyclic amines (harmane, norharmane); "Mel" = melanoidins (A470 / 282 L mol^-1 cm^-1); GO glyoxal; MGO methylglyoxal; AP Amadori product (not measured) |
| Lineage | model system "as previously described (Nguyen et al., 2016; Yu et al., 2020)"; CML/CEL method Jiao 2017 / Yu 2018; dicarbonyls Scheijen & Schalkwijk 2014 / Jiao 2019; fitting Origin 9.0 + SAS 8.0 non-linear least squares (not Athena; no HPD) |

## 1. Why it matters

B18's outcome (`results/validation/kinetic_core_b18_prereg.md` section 6) says the trunk makes almost no
glyoxal from a sugar + amine pot in water at 70-120 C, that its glyoxal sink is a 180 C dry-glass constant
with the barrier fixed at zero, and asks for "the small dicarbonyls in water: their formation from a sugar +
amine pot ... and their loss". This paper measures glyoxal and methylglyoxal in four amine + glucose pots
in 0.1 M phosphate pH 7.0 at both 100 and 130 C over 3-21 min, with lysine, glucose, CML and CEL alongside.
That is exactly the measurement wanted — but the time courses are in Fig. 3 (image) and the text prints
only the glyoxal range for the Lys/Glu pot at each temperature and no methylglyoxal number at all. For
Programme 7 it also gives free-lysine CML and CEL levels at 100 and 130 C in water (Fig. 2D/E, text maxima)
against which a protein-bound estimate can be bounded. Its eleven fitted constants remain unusable: no
units, no rate laws, no intervals, and mislabelled continuation headers (section 3).

## 2. Methods as they matter to a model

- **Buffer:** 0.1 M phosphate, pH 7.0. pH during heating not reported.
- **Charges** (2.2): "Lys–glucose (Lys/Glu; 30 mmol/100 mmol), Asn–glucose (Asn/Glu; 200 mmol/100 mmol),
  Trp–glucose (Trp/Glu; 5 mmol/100 mmol), or a combination of Lys–Asn–Trp–glucose (Mix/Glu; 100 mmol/200
  mmol/5 mmol/100 mmol) were prepared in phosphate buffer". The volume basis is not stated; the modelling
  section says data "were normalized and expressed in millimoles per liter" and the parent system (Nguyen
  2016) is per litre, so read as **mmol/L**: Lys/Glu = 30 mM lysine + 100 mM glucose; Asn/Glu = 200 mM
  asparagine + 100 mM glucose; Trp/Glu = 5 mM tryptophan + 100 mM glucose; **Mix/Glu = 100 mM lysine +
  200 mM asparagine + 5 mM tryptophan + 100 mM glucose**. Note the Mix pot holds 3.3x the lysine of the
  Lys/Glu pot — the paper's "competition" reading of Mix vs Lys/Glu is confounded by this (flag 4). Chosen
  "based on their actual proportions in real cereal and potato-based foods ... appropriately magnified".
- **Vessel / heating:** 10 mL in glass reaction vials with Teflon caps (sealed; vial volume not stated);
  oil bath with magnetic stirrer at **100 C or 130 C**; **5 min pre-heat** to reach temperature, then
  time zero; bath temperature checked every minute, "fluctuations within 3 C"; samples at **3, 6, 9, 12,
  15, 18, 21 min**; ice-quenched; −20 C. Triplicate experiments.
- **Melanoidins:** A470, path 1 cm, ε = 282 L mol^-1 cm^-1 (Zhang 2015) -> "mmol/L" of melanoidin. A
  unit of convenience; not a species concentration.
- **CML and CEL (2.4):** 150 µL sample + 150 µL d4-CML / d4-CEL internal standards; centrifuge, 0.22 µm;
  Waters 2695 + Quattro micro, X-Bridge C18, acetonitrile / 5 mM nonafluoropentanoic acid, ESI+, MRM
  205 -> 84 (CML), 209 -> 88 (d4-CML), 219 -> 84 (CEL), 223 -> 88 (d4-CEL). **Stable-isotope dilution, no
  hydrolysis**: free CML and CEL in solution (there is no protein), directly injected. Internal-standard
  concentrations not printed.
- **GO and MGO (2.5):** 200 µL supernatant + 100 µL 5 µM o-phenylenediamine + 100 µL 4 µM
  2,3-hexanedione internal standard, 4 C dark 12 h; LC-MS/MS MRM of the quinoxalines 131 -> 77 (GO),
  145 -> 77 (MGO), 187 -> 77 (IS). **Not isotope dilution** (non-isotopic dicarbonyl IS). The OPD charge
  (0.5 nmol in 400 µL ≈ 1.25 µM) is far below the measured GO (up to 0.6 mM, section 3) — a derivatising
  agent in this deficit cannot trap the analyte quantitatively unless the sample was diluted first, which
  the text does not say (flag 6).
- **Amino acids, glucose, acrylamide, harmane, norharmane (2.6):** 150 µL sample + 150 µL internal
  standards (200 µg/mL 13C3-acrylamide, 500 µg/mL d4-lysine, 500 µg/mL 13C6-glucose); UHPLC-MS/MS,
  Atlantis dC18, 0.1 % formic acid / acetonitrile, MRM. Lysine, glucose, acrylamide by isotope dilution;
  Asn, Trp, harmane, norharmane against those standards or external.
- **Acrolein:** DNPH derivatisation, HPLC-UV 365 nm. **Acetaldehyde:** headspace GC-FID (75 C 30 min
  equilibration). Both external calibration.
- **Validation (2.9, 3.6, Supplementary Tables 1-2 not on disk):** linearity r > 0.991 over eight
  concentrations; LOD / LOQ at S/N 3 / 10 (values in the supplement); recoveries **87.1-116 % (low
  spike), 83.1-123 % (high spike)**; RSD < 15.3 % intraday, < 14.6 % interday. Which analyte sits where in
  those ranges is not printed.
- **Modelling (2.10, 3.7):** apparent scheme (Fig. 5) keeping reactants (Lys, Trp, Asn, glucose),
  intermediates (acrolein, acetaldehyde, GO, MGO) and products (CML, CEL, harmane, norharmane, acrylamide,
  melanoidins); AP "could not be quantified analytically" and is not a state. "The average duplicate
  analysis data results were normalized and expressed in mmol/L"; k1-k11 by non-linear least squares in
  Origin 9.0 / SAS 8.0 with R^2 per constant; **"the experiments were performed at only one or two
  temperatures; hence, specific activation energies cannot be estimated."** Rate equations: Supplementary
  Fig. 1, absent.
- **Unit conversions used below:** GO 58.04, MGO 72.06, CML 204.22, CEL 218.25, acetaldehyde 44.05,
  acrylamide 71.08 g/mol; 1 µg/mL = 1 mg/L.

### The scheme in words (Fig. 5 is an image; steps inferred from the text's own sentences)

The text names each constant's product; the assignment is quoted so that a later reader can check it:
- **k1 — glyoxal formation** and **k2 — methylglyoxal formation** ("the lower k values for GO, MGO, CML,
  and CEL (k1, k2, k3, and k4)"; "α-dicarbonyl compounds (k1 and k2)"). Source: glucose ("mainly produced
  by thermal glucose degradation and Schiff base decarboxylation") — k1 and k2 are fitted in the Trp and
  Asn pots too, so they are glucose-side steps.
- **k3 — CML formation** (from GO + Lys), **k4 — CEL formation** (from MGO + Lys) — the text's reading of
  "Asn and Trp competing with Lys for the GO and MGO".
- **k5 — acrylamide formation from asparagine** (the Asn route; "CML, CEL, and acrylamide (k3, k4, and
  k5)"; k5 present only in the Asn and Mix pots at 130 C).
- **k6 — acetaldehyde formation** ("acetaldehyde formation (k6)").
- **k7 — acrolein formation from acetaldehyde** ("acetaldehyde oxidation"; k7 is the constant paired with
  k6 in "the k values of acetaldehyde and acrolein").
- **k8 — acrylamide from acrolein** ("the k value of acrylamide (k8) produced by the acrolein oxidation
  pathway"; present in Asn and Mix pots at 130 C only).
- **k9 — harmane**, **k10 — norharmane** ("norharmane (k10)"; k9 and k10 present only where Trp is present
  at 130 C; harmane via Pictet-Spengler condensation of Trp with acetaldehyde, then THβC oxidation).
- **k11 — melanoidins** ("the k11 value of the Lys/Glu model is lower than that of the Mix/Glu model";
  "Mel formation in all models showed the highest k value").
At 100 C "acrylamide, harmane, and norharmane were not detected; hence, the values of k5, k8, k9, and
k10 were considered zero" — matching the dash pattern of Table 1. Orders unknown for every step.

## 3. Tables re-typed

### Table 1. "Estimation of kinetic parameters among different amino acid model systems." (pre-proof pp. 32-33)

Printed header, first block: "K1 (10–3) R2 | K2 (10–2) R2 | K3 (103) R2 | K4 (102) R2 | K5 (102) R2".
Printed header, "Table 1. Continued": "K1 (10–3) R2 | K2 (10–2) R2 | K3 (103) R2 | K4 (102) R2 | K5
(102) R2 | K1 (10–3) R2" — **a copy of the first block's header with a sixth column added; the six
continuation columns are k6-k11** (the caption says k1-k11; the dash pattern matches the text's k5/k8/k9/k10
absences and the k6/k7/k10/k11 comparisons quoted in section 2). **The scale prefixes of k6-k11 are
therefore unknown**, and the prefixes of k1-k5 are ambiguous in direction (does "K3 (10^3)" mean the
printed number is k3 x 10^3, i.e. k3 = 2.94 x 10^-3, or that k3 is in thousands?). **No units are printed
for any constant.** Footnote a: "Data were expressed as mean ± SD in triplicates (n = 3) and fit by the
non-linear curve fitting method" — **no ± SD appears in any cell.**

| pot / T | k1 (GO) | R^2 | k2 (MGO) | R^2 | k3 (CML) | R^2 | k4 (CEL) | R^2 | k5 (acrylamide, Asn) | R^2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Lys/Glu 100 C | 8.34 | 0.91 | 0.66 | 0.95 | 2.94 | 0.82 | 0.25 | 0.81 | – | |
| Trp/Glu 100 C | 0.02 | 0.78 | 0.20 | 0.96 | – | | – | | – | |
| Asn/Glu 100 C | 1.80 | 0.93 | 0.16 | 0.99 | – | | – | | – | |
| Mix/Glu 100 C | 5.58 | 0.89 | 0.68 | 0.95 | 1.19 | 0.90 | 0.19 | 0.93 | – | |
| Lys/Glu 130 C | 8.15 | 0.76 | 3.70 | 0.93 | 6.11 | 0.81 | 1.25 | 0.91 | – | |
| Trp/Glu 130 C | 3.15 | 0.91 | 1.55 | 0.94 | – | | – | | – | |
| Asn/Glu 130 C | 4.43 | 0.79 | 3.16 | 0.95 | – | | – | | 6.05 | 0.81 |
| Mix/Glu 130 C | 4.67 | 0.83 | 4.05 | 0.89 | 1.91 | 0.83 | 0.77 | 0.93 | 1.09 | 0.89 |

| pot / T | k6 (acetaldehyde) | R^2 | k7 (acrolein) | R^2 | k8 (acrylamide via acrolein) | R^2 | k9 (harmane) | R^2 | k10 (norharmane) | R^2 | k11 (melanoidins) | R^2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Lys/Glu 100 C | 0.27 | 0.88 | 1.50 | 0.77 | – | | – | | – | | 1.41 | 0.98 |
| Trp/Glu 100 C | 0.16 | 0.77 | 0.94 | 0.79 | – | | – | | – | | 0.11 | 0.83 |
| Asn/Glu 100 C | 0.15 | 0.81 | 0.63 | 0.75 | – | | – | | – | | 0.59 | 0.99 |
| Mix/Glu 100 C | 0.19 | 0.80 | 0.79 | 0.79 | – | | – | | – | | 1.76 | 0.98 |
| Lys/Glu 130 C | 1.58 | 0.98 | 11.4 | 0.79 | – | | – | | – | | 1.58 | 0.94 |
| Trp/Glu 130 C | 1.09 | 0.96 | 9.58 | 0.80 | – | | 4.47 | 0.78 | 8.35 | 0.87 | 0.65 | 0.97 |
| Asn/Glu 130 C | 5.77 | 0.98 | 445 | 0.77 | 4.41 | 0.92 | – | | – | | 1.21 | 0.95 |
| Mix/Glu 130 C | 5.03 | 0.92 | 405 | 0.83 | 4.02 | 0.88 | 4.79 | 0.77 | 7.84 | 0.90 | 1.66 | 0.92 |

Checks against the text (all with my column assignment): "k1 ... Mix lower than Lys/Glu" 5.58 < 8.34
(100 C), 4.67 < 8.15 (130 C) — yes; "k2 lower in Mix" 0.68 vs 0.66 and 4.05 vs 3.70 — **no, higher at
both temperatures** (the text's claim is contradicted by its table); k3, k4 lower in Mix — yes; k6 Mix >
Trp at 130 C (5.03 vs 1.09) — yes; k6, k7 Mix < Asn (5.03 vs 5.77; 405 vs 445) — yes; k8 Mix < Asn (4.02
vs 4.41), k10 Mix < Trp (7.84 vs 8.35) — yes; k11 Lys/Glu < Mix at 100 C (1.41 vs 1.76) — yes. "R^2
ranging from 0.80-0.99" — the table has 0.75-0.79 in eleven cells. "k values generally increased with
temperature" — k1 in Lys/Glu falls (8.34 -> 8.15). R^2 columns are per-constant fits, which means the
constants were fitted one response at a time, not as a joint multiresponse estimate.

### Concentrations printed in the text (Figs 1-4 are figure-only)

| quantity | value | mmol/L (mine) | condition | source |
|---|---|---:|---|---|
| Lys lost by 21 min | 29.1 % (Lys/Glu), 19.3 % (Mix) | 8.7 of 30; 19 of 100 | 100 C | 3.1 |
| Lys lost by 21 min | 56.1 % (Lys/Glu), 38.7 % (Mix) | 16.8 of 30; 39 of 100 | 130 C | 3.1 |
| Trp lost | 20.7 % / 6.0 % (Trp/Glu / Mix) at 100 C; 79.2 % / 10.9 % at 130 C | | | 3.1 |
| Asn lost | 33.3 % / 30.9 % at 100 C; 52.0 % / 50.1 % at 130 C | | | 3.1 |
| glucose consumed at 100 C | 15.7 % (Lys/Glu), 28.1 % (Mix), 24.1 % (Asn/Glu), 15.8 % (Trp/Glu) | 15.7 of 100 mM | 100 C, 21 min | 3.1 |
| glucose consumed at 130 C | 63.1 % (Lys/Glu), 95.5 % (Mix), 73.2 % (Asn/Glu), 66.1 % (Trp/Glu) | 63 of 100 mM | 130 C, 21 min | 3.1 |
| **GO, Lys/Glu, 100 C** | **3.01 ± 0.50 to 7.38 ± 1.54 µg/mL** over the run | **0.052 to 0.127** | 100 mM Glc + 30 mM Lys, pH 7.0 | 3.4 |
| **GO, Lys/Glu, 130 C** | **8.33 ± 1.81 to 35.1 ± 0.43 µg/mL** | **0.144 to 0.605** | same | 3.4 |
| GO shape | rises with time in all pots except Trp/Glu at 100 C and Lys/Glu at 130 C, which "showed a maximum GO concentration at the initial stage, but significantly decreased"; Lys/Glu highest, Mix next | | | 3.4 |
| MGO | "significant increase ... with increased reaction time, reaching a relatively high value at the end"; Mix and Lys/Glu highest at 100 C; Mix and Asn/Glu highest at 130 C; **no number printed** | — | | 3.4, Fig. 3B/D |
| CML, Lys/Glu | 11.6 ± 2.21 µg/mL (100 C); 27.2 ± 2.31 µg/mL (130 C) | 0.057; 0.133 | maxima (Fig. 2D) | 3.3 |
| CML, Mix | 5.11 ± 1.03 (100 C); 9.06 ± 1.89 µg/mL (130 C) | 0.025; 0.044 | | 3.3 |
| CEL, Lys/Glu | 10.6 ± 0.49 µg/mL (100 C); 51.7 ± 8.74 µg/mL (130 C) | 0.049; 0.237 | | 3.3 |
| CEL, Mix | 8.30 ± 3.15 (100 C); 31.9 ± 8.09 µg/mL (130 C) | 0.038; 0.146 | | 3.3 |
| acetaldehyde, Lys/Glu | 3.02 ± 0.12 (100 C); 13.1 ± 1.13 µg/mL (130 C); Asn/Glu 12.2 ± 1.53 at 130 C | 0.069; 0.297 | | 3.5 |
| acrolein at 130 C | Lys/Glu 0.06 ± 0.01 to 1.29 ± 0.12; Trp/Glu 0.20 ± 0.09 to 1.25 ± 0.09 µg/mL | 0.001-0.023 | | 3.5 |
| acrylamide at 130 C | Asn/Glu 22.9 ± 2.55; Mix 12.6 ± 2.80 µg/mL (max); none at 100 C; appears after 6 min | 0.322; 0.177 | 200 mM Asn + 100 mM Glc | 3.2 |
| norharmane at 130 C | Trp/Glu 27.9 ± 2.03; Mix 22.6 ± 1.74 ng/mL (max, from 6 min) | 1.7e-4; 1.3e-4 | | 3.2 |
| harmane at 130 C | Trp/Glu 13.5-13.9 ng/mL (from 15 min); Mix +6.6-13.1 % | 8e-5 | | 3.2 |
| melanoidins, Mix | 3.37 ± 0.94 (100 C); 30.8 ± 0.75 "mmol/L" (130 C); Trp/Glu lowest | (A470-based) | 21 min | 3.3 |

Molar reading of the AGE data in the Lys/Glu pot (mine): at 100 C, 21 min, lysine lost 8.7 mM; CML 0.057 +
CEL 0.049 = 0.106 mM = **1.2 % of the lysine lost**; glyoxal present 0.05-0.13 mM (~0.1 % of glucose). At
130 C: lysine lost 16.8 mM; CML + CEL = 0.37 mM = **2.2 % of the lysine lost**; glyoxal up to 0.6 mM
(0.6 % of glucose). CEL ≥ CML at both temperatures (0.9 : 1 at 100 C, 1.8 : 1 at 130 C) — the same
direction as Berk 2021 in sesame (2.9 : 1) and the opposite of Nguyen 2016's casein (CML > CEL).

## 4. Kinetic numbers the repository can use

Registry keys: CML -> `cml`; CEL -> `cel`; acrylamide -> `acrylamide`; acetaldehyde -> `acetaldehyde`;
acrolein -> `acrolein`; glyoxal, methylglyoxal, lysine, harmane, norharmane: **not in registry**. Conditions for every row: 0.1 M phosphate pH 7.0, sealed stirred vial, 10 mL, 5-min
pre-heat then 3-21 min, triplicate.

| step | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| glyoxal in a glucose + lysine pot | level range over 3-21 min | 0.052-0.127 | mmol/L | 100 C, 100 mM Glc + 30 mM Lys | 3.4 | level_only |
| glyoxal in a glucose + lysine pot | level range over 3-21 min | 0.144-0.605 | mmol/L | 130 C, same | 3.4 | level_only |
| glyoxal, 130 / 100 C | ratio of upper values | 4.8 | — | Lys/Glu | derived from 3.4 | within_study_ratio |
| glyoxal shape, Lys/Glu 130 C | early maximum then decline | — | | 3.4 | structural (loss is visible at 130 C in water, formation dominates at 100 C) |
| glyoxal and methylglyoxal time courses, four pots, two temperatures | — | mmol/L vs min | | Fig. 3 | figure_only |
| methylglyoxal | rises monotonically to 21 min in all pots | — | 100 and 130 C | 3.4 | structural, no number |
| CML from free lysine | maximum | 0.057 / 0.133 | mmol/L | 100 / 130 C, Lys/Glu | 3.3 | level_only |
| CEL from free lysine | maximum | 0.049 / 0.237 | mmol/L | 100 / 130 C, Lys/Glu | 3.3 | level_only |
| (CML + CEL) / lysine lost | 0.012 / 0.022 | — | 100 / 130 C, Lys/Glu, 21 min | derived | within_study_ratio |
| CEL / CML | 0.9 / 1.8 | — | 100 / 130 C, Lys/Glu | derived | within_study_ratio |
| CML: Mix / Lys-Glu | 0.44 / 0.33 | — | 100 / 130 C (Mix has 3.3x the lysine — flag 4) | derived | within_study_ratio (confounded) |
| lysine loss, Lys/Glu | 29.1 % / 56.1 % in 21 min | — | 100 / 130 C | 3.1 | level_only; apparent first-order 1.6e-2 / 3.9e-2 min^-1 (mine, single-point) |
| glucose loss, Lys/Glu | 15.7 % / 63.1 % in 21 min | — | 100 / 130 C | 3.1 | level_only |
| acetaldehyde | 0.069 / 0.297 | mmol/L | 100 / 130 C, Lys/Glu | 3.5 | level_only |
| acrylamide from Asn + Glc | 0.322 max | mmol/L | 130 C, 200 mM Asn + 100 mM Glc; none at 100 C | 3.2 | level_only (already in the acrylamide family's purview) |
| k1-k11 (Table 1) | see section 3 | **no units, no orders, no intervals** | | Table 1 | **refused** (k3 C.5 stands) |

What B18 gets from this paper: an order of magnitude for glyoxal in water from a sugar + amine pot — about
0.1 mM at 100 C and 0.6 mM at 130 C from 100 mM glucose within 20 min, with visible glyoxal LOSS only at
130 C in the lysine pot. The B18 note that the trunk "makes almost no glyoxal from a sugar + amine pot" can
be tested against these two levels as a level-only check (a factor-of-several tolerance, given flag 6). No
rate constant for glyoxal formation or loss can be taken.

## 5. Flags

1. **The rate-constant table is unusable, as already declared** (k1 §2d, k3 C.5, FIT_HOLDOUT_DECLARATION):
   no units anywhere; the prefixes' direction is ambiguous; the continuation header is a copy of the first
   block so k6-k11 have no printed prefix at all; no SD despite the footnote; rate equations in an absent
   supplement, so no reaction orders; per-constant R^2 means single-response fits, not multiresponse. The
   column-to-step mapping in section 3 is mine, built from the text's own sentences, and it is consistent
   with every quoted comparison except the k2 claim; it does not make the numbers usable.
2. **Pre-proof.** The version of record may have corrected the Table 1 headers; it is not on disk (also
   asked for in `k3_final_parameter_inventory.md` row 10 of the wishlist).
3. **All time courses are figure-only.** For glyoxal only the Lys/Glu range at each temperature is
   printed; for methylglyoxal nothing. The maxima/ranges quoted are what the repository can hold.
4. **The Mix/Glu pot is not the Lys/Glu pot plus Asn and Trp**: it has 100 mM lysine against 30 mM. Every
   "competition" ratio (CML, CEL, k3, k4) between Mix and Lys/Glu is confounded by a 3.3x lysine change and
   by the pH-7 buffer's capacity against 305 mM of amino acids. Use Mix numbers as a separate pot.
5. **Time base**: t = 0 is after a 5-min pre-heat at bath temperature; the first sample is at 3 min (8 min
   of heating). Any comparison with a model that starts cold must add the pre-heat.
6. **Dicarbonyl derivatisation stoichiometry**: 100 µL of 5 µM OPD into 200 µL sample (+ 100 µL internal
   standard) gives ~1.25 µM OPD against up to 600 µM glyoxal plus methylglyoxal; unless the sample was diluted (not stated) the
   quinoxaline yield is reagent-limited and the GO/MGO levels are lower bounds with a non-isotopic internal
   standard. The Scheijen & Schalkwijk 2014 source method uses OPD in large excess; the printed "5 µM" may
   be a slip for 5 mM. Treat the absolute GO levels as ±0.5 dex until the version of record is checked.
7. **Melanoidin "mmol/L"** is A470 / 282 — an absorbance in molar disguise; not a concentration to balance.
8. **Recoveries 83-123 %, RSD to 15 %**; which analyte is at which end is in the absent supplement.
9. **No pH drift, no headspace volume, no oxygen statement**; sealed Teflon-capped vials, stirred.
10. **Free lysine, not protein-bound**: for Programme 7 these CML/CEL levels are an upper-bound reference
    (free ε-amines at 30 mM), not a protein datum; Nguyen 2016 (casein) and Troise 2015 (soy) are the
    protein-bound comparators.

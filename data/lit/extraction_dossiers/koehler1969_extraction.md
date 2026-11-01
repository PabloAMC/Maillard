# Koehler, Mason & Newell 1969 — EXTRACTION
DOI: 10.1021/jf60162a006 (printed: none; Crossref 2026-11-01 resolves to "Formation of pyrazine compounds in sugar-amino acid model systems", J. Agric. Food Chem. 17, 393-396, 1969, authors Koehler / Mason / Newell)
### *"Formation of Pyrazine Compounds in Sugar-Amino Acid Model Systems"*
### Glucose (or sucrose / fructose) + asparagine, glutamine, aspartic acid, glutamic acid, alanine or NH4Cl, 10 mmol each in 20 mL water + 200 mL diethylene glycol, reflux 120 °C, 24 h; falling-film isolation, packed-column GC. Pyrazine distributions per N source (Fig. 1, figure-only), µmol per product with and without excess acetaldehyde (Fig. 3, figure-only), total yields (Table I), and 14C origin of MPZ / DMP carbon (Tables II-III). **No pH measured; "dimethylpyrazine" on a packed Carbowax column (2,5 / 2,6 not shown resolved).**

**Read 2026-11-01 from `data/articles/koehler1969.pdf`** (4 pp., JAFC 17(2):393-396; Acrobat Capture scan with an OCR text
layer). Prose from `pdftotext -layout`; Tables I-III checked against a 200 dpi render; **Figs. 1 and 3 read off 400 dpi
crops** (bar tops; read error ±1 percentage point on Fig. 1, ±0.3 point on bars under 5 %; ±3-8 µmol on Fig. 3). Every
derived number was computed in the `maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | P. E. Koehler, M. E. Mason, J. A. Newell (Dept. of Biochemistry, Oklahoma State University, Stillwater) |
| Venue | J. Agric. Food Chem. 17 (2), 393-396, Mar-Apr 1969; received 26 Jun 1968 |
| Laboratory | Oklahoma State (Koehler / Mason / Newell / Odell). **Not in `docs/reference/contrast_benchmark/lab_ids.yml`**: a new lab id is needed if any cell is keyed. Independent of UC Davis, Reineccius and the Ho network |
| Companion | koehler1970 (same lab, same solvent and isolation, 10x the reactant concentration) |

## 1. What this paper answers for the model

Context: EX-052 (`docs/history/pyrazine_rivals_2026-10-31/`). Mixing ratio R = MPZ² / (PZ x DMP): one shared constant
with the statistical factor (P1) caps R at 4 (x h, the pool-heterogeneity factor <= 1); the engine as coded gives R = 1
(CF-077); pair-specific constants (P1x) or pair-specific dihydropyrazine fate (P1y) allow any R >= 4h; an extra MPZ
route (P2b) allows any R.

| need | what Koehler 1969 gives | tag |
|---|---|---|
| PZ / MPZ / DMP per pot | **Fig. 1** (% of total pyrazines, glucose + Asn / Gln / Glu / Asp / NH4Cl, and roasted peanuts); **Fig. 3** (µmol, glucose + Asn ± excess acetaldehyde). Figure-only | [F] |
| R | Asn 23 (Fig. 1) or **50** (Fig. 3, band 34-82); Gln 39; Asp 14; Glu 9.3; **NH4Cl 0.07** (band 0.015-0.29); roasted peanuts ~46 (band 23-300); Asn + acetaldehyde **98** (band 72-145) | [Z] |
| response factors | "peak areas of samples with those of known weights of standards" [V]: per-compound external calibration on FID. Recovery > 90 % [M]. Whether Fig. 1's "percent yield" is a mole or weight share is not stated (R differs by 2 % only) | [M] / [!] |
| 2,5- vs 2,6-DMP | not resolved or not reported: packed 20 ft Carbowax 20M; Fig. 1 draws one "2,5-dimethyl" bar. Treat "DMP" as 2,5(6) | [!] |
| pH dependence | none (glycol-water; pH not measured) | — |
| temperature dependence | none (120 °C only) | — |
| nitrogen source | **yes, the paper's main within-study factor**: four amino acids vs NH4Cl give different distributions (Fig. 1) and total yields (Table I) | [M] / [F] |
| 13C / 14C carbon origin | MPZ / DMP carbon from glucose; Asn-UL-14C gives MPZ / DMP specific-activity ratio 0.002 / 0.007, Ala-2-14C 0.026 / 0.007 (Table II) | [M] |
| separates P1x / P1y from P2b? | **No** (sec. 5). The acetaldehyde arm is a supply perturbation under which pairing predicts R unchanged; R rose ~2x (50 -> 98), but the read-error bands nearly touch | [Z] |

**Answers.**

1. **R >> 4 for every amino-acid pot (9-50; a fourth laboratory beside Reineccius, the Ho network and UC Davis),
   R ~0.07 for NH4Cl.** The amino-acid values sit in Leahy's range (7-70), far below glucose + NH3 at UC Davis
   (300-700). The NH4Cl pot is the opposite extreme: pyrazine is ~75 % of the total.
2. **Amino-acid carbon in MPZ <= 2.6 % (Ala-2-14C) and <= 0.2 % (Asn-UL-14C)** [M, Table II]. That is a third
   laboratory (with Nestlé and Reading, S4) against an aldehyde-addition MPZ route that carries amino-acid carbon
   (P2a_aa), and at 120 °C in glycol-water rather than 180 °C roasts. It says nothing about a sugar-carbon extra route (P2b).
3. **Glucose C1 in MPZ: 1.29 / 1.20 (Asn / Ala) vs 1.50 for "one C2 (C1-C2) + one C3" pairing; C6 0.63 and C3,4
   0.63 vs 0.50; DMP 0.97-1.14 vs 1.00** [M, Table III]. The authors read the shortfall as dilution of the C2 pool by
   further C4 cleavage [Q]. The same numbers fit an MPZ fraction formed from a non-C1-C2 C2 unit or by a route other than
   C2 + C3 pairing. Not separable here (the C2-pool composition is free).

## 2. Conditions [M]

| item | value |
|---|---|
| Reactants | 10 mmol amino acid (or NH4Cl) + 10 mmol sugar (glucose; also sucrose, fructose in Table I) |
| Medium | 20 mL deionised water + 200 mL diethylene glycol ("low water content"); ~45 mM each reactant [Z] |
| pH / buffer | none; not measured [!] |
| T / time | reflux at 120 °C, 24 h, stirred; open reflux (air) |
| Isolation | falling-film evaporator at 98 °C, liquid-N2 cold finger, CH2Cl2 extraction; recovery "better than 90%" [V] |
| GC | Perkin-Elmer 801 FID, 20 ft x 1/4 in glass, 15 % Carbowax 20M, 75-190 °C at 4 °C/min; identity by retention vs standards, MS (LKB 9000 prototype) checked on the Asn system |
| Total yield | UV absorbance at λmax ~272-280 nm calibrated with 2,5-DMP (ε 6.97 absorbance units per µmol [M]); a lumped estimate |
| Labelling | 2 µCi of 14C-glucose (C1, C6, C3,4) or 14C-amino acid; MPZ and DMP trapped by preparative GC (TCD), scintillation counting; Table III values are means of 2-3 runs ± average deviation |
| Replicates | Figs. 1, 3: not stated; no error bars |

## 3. Quantitative results

**Table I (total pyrazine yield, µmol, UV estimate) [M]:** Asn-glucose 411; Gln-glucose 205; Asp-glucose 198; Glu-glucose
114; NH4Cl-glucose 59; NH4Cl-fructose 195 (300 dpi render); Asn-sucrose 268; Asn-glucose-acetaldehyde **11643**
(render confirms the digits) [!]: 28x the Asn pot, while Fig. 3's bars sum to ~580 µmol. The UV estimate is
calibrated on 2,5-DMP and acetaldehyde condensates absorb near 280 nm, so the 11643 is probably not pyrazines alone [Z].
λmax 262 (NH4Cl-glucose) to 280 nm. Glucose-Asn ≈ 4 % of glucose [M].

**Fig. 1 (% of total pyrazines; bar tops read at 400 dpi) [F]:**

| N source (+ glucose) | PZ | MPZ | DMP (2,5(6)) | R [Z] |
|---|---:|---:|---:|---:|
| roasted peanuts (distillate) | ~0.6 | ~29.7 | ~31.6 | ~46 (23-300) |
| Asn | ~4.8 | ~54.9 | ~26.6 | 23 (20-28) |
| Gln | ~4.8 | ~61.3 | ~19.8 | 39 (34-47) |
| Glu | ~14.4 | ~31.6 | ~7.3 | 9.3 (8-11) |
| Asp | ~2.8 | ~28.4 | ~20.9 | 13.5 (11-18) |
| NH4Cl | ~74.5 | ~2.3 | ~1.0 | 0.07 (0.015-0.29) |

Glu gives ~40 % of a trialkyl bar (2-ethyl-3,6(5)-dimethyl, by the drawn structure) [F].

**Fig. 3 (µmol recovered, glucose + Asn ± "a large excess" of acetaldehyde, amount not stated [!]) [F]:**

| pot | PZ | MPZ | DMP | last bar (TrMP-type) | R [Z] |
|---|---:|---:|---:|---:|---:|
| Asn | ~9 (±3) | ~224 | ~112 | ~17 | 50 (34-82) |
| Asn + acetaldehyde | ~157 | ~401 | ~10.5 | ~10 | 98 (72-145) |

Acetaldehyde: PZ x~17, MPZ x~1.8, DMP x~0.09 [Z]. Fig. 1's Asn PZ share (4.8 %) and Fig. 3's (9 / 362 = 2.5 %) disagree
by ~2x inside the paper [Z][!]; R moves with it (23 vs 50).

**Tables II-III (ratio of product to reactant specific activity) [M]:**

| labelled reactant | partner | MPZ | DMP |
|---|---|---:|---:|
| glucose-1-14C | Asn | 1.290 | 0.970 |
| Asn-UL-14C | glucose | 0.002 | 0.007 |
| glucose-1-14C | Ala | 1.200 | 0.950 |
| Ala-2-14C | glucose | 0.026 | 0.007 |
| glucose-1-14C (Table III) | Asn | 1.29 ± 0.01 (pathway I 1.50, II 0.83) | 0.97 ± 0.01 |
| glucose-6-14C | Asn | 0.63 ± 0.01 (I 0.50, II 0.83) | 1.14 ± 0.20 |
| glucose-3,4-14C | Asn | 0.63 ± 0.005 (I 0.50, II 0.83) | 1.10 |

## 4. Mechanism evidence and its kind

- **Labelling (kind 1, 14C), one lab (Oklahoma State):** sugar carbon in MPZ and DMP; amino-acid carbon <= 2.6 %.
  Position data favour a single C2 unit from glucose C1-C2 and C3 units from both halves [Q].
- **Within-study factor (kind 6):** N source (Fig. 1), acetaldehyde addition (Fig. 3).
- The C2 + C3 fragment pairing, ammonolysis of a C6 intermediate, and the claim that ammonia is not the common N carrier
  are **interpretations** [Q]/[P]; the NH4Cl contrast is confounded by pH (NH4Cl is acidic, amino acids near neutral).

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 P0 (engine, R = 1) / CF-077 | CONTRADICTS P0 (fourth lab); CF-077 unchanged | amino-acid pots R 9-50 |
| EX-052 P1 (R <= 4h) | CONTRADICTS (amino-acid pots) | R/4 = 2-12; figure-read, no RF check printed per compound, 2,5 / 2,6 unresolved (co-elution lowers R, so the violation is conservative) |
| P1x / P1y vs P2b | **LEAVES UNSEPARATED** | Under any pairing rival with fixed constants R does not depend on the A : M supply; adding acetaldehyde (a C2 / A-side perturbation) moved R 50 -> 98 (x2; bands 34-82 vs 72-145). Direction fits P2b or a heterogeneity change, magnitude within ~1.5 read errors. Weak, single figure, amount of acetaldehyde unknown |
| S4 / P2a_aa (amino-acid-carbon MPZ) | SUPPORTS rejection (third lab) | Ala-2-14C: 2.6 % in MPZ; Asn-UL: 0.2 % |
| low-R rows (CF-087, the Ho glycine row) | bears on it | NH4Cl-glucose R ~0.07: a second lab and system with R < 0.4h_any, i.e. a row every pairing rival forbids unless h < 0.02 or extra PZ / DMP routes exist |
| evidence map 2.7 "methylpyrazine extra route: NONE" | LEAVES UNCHANGED | no labelled discrimination of an extra sugar-carbon MPZ route |
| E4 (pyrazine pH law) | LEAVES UNCHANGED | no pH series |

## 6. Disagreements found

1. **R in glucose + inorganic N, two labs** [Z]: NH4Cl-glucose here R ~0.07 (PZ ~75 % of pyrazines) vs glucose + NH4OH at
   UC Davis R 310-700 (MPZ ~80 %; shibamoto1976 / 1977). ~3.7 dex. Candidate causes: pH (NH4Cl acidic vs 8 M NH3
   strongly basic; Koehler 1970 Fig. 4 shows acid suppresses pyrazines), medium (glycol-water vs water), T (120 vs
   100 °C). Koehler 1970 itself found NH4Cl "ND" at 10x concentration (disagreement 2).
2. **Within lab, NH4Cl-glucose: 59 µmol here vs not detected in koehler1970 Table I** (10x concentration, same
   solvent and T). Cause unknown (the 1970 paper notes NH4Cl results "similar" to aspartate's acidity effect).
3. **Within lab, Asn-glucose distribution: DMP/MPZ ~0.5 here (Fig. 1 0.48, Fig. 3 0.50) vs 3.7 in koehler1970
   Table I** (same 120 °C / 24 h / solvent; 10 vs 100 mmol in 220 mL). Candidate causes: reactant concentration (supply
   ratio and second-order steps), the 1970 isolation (two falling-film passes), the DMP peak's identity on the packed
   column. Koehler 1970 Fig. 2 has MPZ > DMP only before ~3 h, so time is not the cause.
4. **Inside the paper: Asn PZ share 4.8 % (Fig. 1) vs ~2.5 % (Fig. 3)**, so R 23 vs 50.

## 7. Ledger-ready measurements

1. Table II (14C specific-activity ratios, MPZ and DMP, Asn and Ala) — the amino-acid-carbon bound for S4 from a
   third lab, printed numbers.
2. Table III (glucose C1 / C6 / C3,4 into MPZ and DMP) — fragment-origin data, printed with deviations.
3. Fig. 3 (µmol, ± acetaldehyde) — the only supply perturbation with all three species; figure-only, two readers.
4. Fig. 1 NH4Cl and Asn bars — R extremes; figure-only, two readers; key DMP as `2,5(6)-DMP unresolved`.

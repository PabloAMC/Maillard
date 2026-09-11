# De Vleeschouwer, Van der Plancken, Van Loey & Hendrickx 2009 — Part I: effect of the type of sugar — EXTRACTION (low-moisture asparagine–sugar acrylamide kinetics, 120–200 °C, multiresponse, ± 95 % HPD)

**Source on disk:** `data/articles/devleeschouwer2009b.pdf` (Food Chemistry 114 (2009) 116–126;
DOI 10.1016/j.foodchem.2008.09.024). Read 2026-09-11 with `pypdf` (59 k characters, full text layer;
Tables 1–3 transcribed below from the text layer; Schemes 1–4 and Figures are images). Dossier
written because the 2026-09-11 audit found this the one PDF on disk that another dossier
(`devleeschouwer2009_extraction.md`, **Part II**, flag 2) explicitly declares "a DIFFERENT PAPER and
it has no dossier".

| field | value |
|---|---|
| Title | "Role of precursors on the kinetics of acrylamide formation and elimination under low moisture conditions using a multiresponse approach – Part I: Effect of the type of sugar" |
| Authors | Kristel De Vleeschouwer, Iesel Van der Plancken, Ann Van Loey, Marc E. Hendrickx (KU Leuven, LFoRCe) |
| Systems | **equimolar asparagine–sugar**, sugar = glucose, fructose or sucrose; equilibrated at initial a_w **0.92 (at 4 °C)**; initial moisture (Table 1): glucose 14.53 ± 0.08 %, fructose 12.60 ± 0.12 %, sucrose 17.79 ± 0.28 % |
| Conditions | 120, 140, 160, 180, 200 °C; time courses; **low moisture**, not aqueous |
| Responses | acrylamide (GC-MS, CI), glucose / fructose / sucrose (HPAEC-PAD), asparagine and aspartic acid (EZ:faast GC-MS), melanoidins (A470, ε = 282 L mol⁻¹ cm⁻¹, Knol 2005) |
| Model | multiresponse Bayesian fit, reparameterised Arrhenius with T_ref = 160 °C; Scheme 4 for glucose and fructose, Scheme 3 for sucrose (adds hydrolysis k_HY and inversion k_I) |
| Companion | `devleeschouwer2009_extraction.md` = **Part II** (effect of the type of amino acid / other precursors), same laboratory and method |

## 1. Table 2 — relative maximum acrylamide yield per mol initial asparagine (%), glucose = 100

| T (°C) | glucose | fructose | sucrose |
|---:|---:|---:|---:|
| 120 | 100 | 165.3 | 220.3 |
| 140 | 100 | 152.9 | 184.3 |
| 160 | 100 | 110.3 | 151.8 |
| 180 | 100 | 148.6 | 163.8 |
| 200 | 100 | 127.8 | 116.0 |

**Within-study ratios** (fructose/glucose, sucrose/glucose at five temperatures) — the class of number
this repository FITS on. Sucrose's excess collapses from 2.2× at 120 °C to 1.16× at 200 °C.

## 2. Table 3 — kinetic parameters, T_ref = 160 °C, ± 95 % HPD (verbatim; "–" = not in the model)

| parameter | glucose (Scheme 4) | fructose (Scheme 4) | sucrose (Scheme 3) |
|---|---:|---:|---:|
| k_F,ref (10⁻³ mM⁻¹ min⁻¹), acrylamide formation | 3.57 ± 1.38 | 7.40 ± 9.48 | 3.56 ± 0.86 |
| k_E,ref (min⁻¹), acrylamide elimination | 0.10 ± 0.04 | 0.09 ± 0.02 | 0.74 ± 0.17 |
| k_INTg,ref (M⁻¹ min⁻¹) | 1.70 ± 1.05 | – | 1.28 ± 1.09 |
| k_INTf,ref (M⁻¹ min⁻¹) | – | 0.22 ± 0.38 | 1.73 ± 0.80 |
| k_M,ref (min⁻¹), melanoidin | 1.23 ± 0.49 | 0.58 ± 0.16 | 0.04 ± 0.01 |
| k_B,ref (M⁻¹ min⁻¹) | 3.90 ± 3.68 | 0.63 ± 30.64 | 4.49 (fixed, footnote c) |
| k_Cg,ref (10⁻³ min⁻¹) | indeterminate | – | 0.48 ± 0.23 |
| k_Cf,ref (10⁻³ min⁻¹) | – | 1.03 ± 0.86 | 716.18 ± 122.10 |
| k_Asp,ref (10⁻³ min⁻¹), aspartic acid | 26.43 ± 5.76 | 13.62 ± 4.08 | 7.29 ± 2.60 |
| k_X,ref (10⁻³ min⁻¹) | indeterminate | 2.12 ± 3.25 | 1.28 ± 5.11 |
| k_HY,ref (10⁻³ min⁻¹), sucrose hydrolysis | – | – | 0.47 ± 0.09 |
| k_I,ref (10⁻³ min⁻¹), inversion | – | – | 0.50 ± 0.20 |
| Ea_F (kJ/mol) | 159.2 ± 29.5 | 122.0 ± 66.5 | 96.2 ± 22.3 |
| Ea_E (kJ/mol) | 113.2 ± 32.3 | 95.4 ± 18.3 | 108.9 ± 21.8 |
| Ea_INTg (kJ/mol) | 117.5 ± 25.2 | – | 328.2 ± 61.3 |
| Ea_INTf (kJ/mol) | – | 149.1 ± 87.7 | 180.7 ± 21.1 |
| Ea_M (kJ/mol) | 105.7 ± 29.1 | 90.1 ± 27.8 | 64.7 ± 17.8 |
| Ea_B (kJ/mol) | 180.3 ± 38.5 | 119.6 ± 44.3 | 34.7 (fixed) |
| Ea_Cg (kJ/mol) | −6.7 ± 0.2 | – | indeterminate |
| Ea_Cf (kJ/mol) | – | 152.7 ± 44.1 | 124.3 ± 10.5 |
| Ea_Asp (kJ/mol) | 105.4 ± 10.6 | 108.3 ± 11.4 | 109.4 ± 16.1 |
| Ea_X (kJ/mol) | 668.9 ± 35.2 | 167.6 ± 68.8 | 322.2 ± 296.0 |
| Ea_HY (kJ/mol) | – | – | 140.6 ± 8.5 |
| Ea_I (kJ/mol) | – | – | 113.2 ± 13.8 |

The paper marks parameters "in italic" as having no physical meaning (not recoverable from the text
layer; Ea_X at 668.9 and 322.2 and Ea_INTg at 328.2 are the obvious candidates). Footnote c: fixed at
their estimated values and the rest re-estimated.

## 3. What this is for this model

- **Lane: acrylamide.** Rates and barriers with HPDs on the acrylamide formation and elimination
  steps from asparagine + glucose at 120–200 °C — FIT-class evidence under the standing rule. The
  acrylamide lane's fit (`kinetic_core_b3_fit_report.json`) does not read this paper; whether it should
  is a **pre-registration decision**, not made here. The a_w is 0.92 initial in a low-moisture matrix
  (14 % moisture), which is inside the lane's declared a_w window but is not an aqueous pot.
- **Cross-check available today, no fit:** the lane's shipped `log10_k_ref_at_160C` and
  `fitted_Ea_kJ_mol` for the Asn + Glc initiation and the elimination step can be compared with
  k_F,ref = 3.57e-3 mM⁻¹ min⁻¹ / Ea_F = 159.2 kJ/mol and k_E,ref = 0.10 min⁻¹ / Ea_E = 113.2 kJ/mol
  (glucose column). Different reference temperature conventions must be reconciled first.
- **Fructose and sucrose columns** are the type-of-sugar directional claims: fructose gives more
  acrylamide than glucose at every temperature (Table 2), a claim the trunk's fructose route can be
  asked about.

## 4. Flags

1. The text layer loses italics, so which parameters the authors call physically meaningless is
   inferred, not read.
2. Schemes 3 and 4 (the reaction networks that define k_INT, k_M, k_B, k_C, k_X) are images; the
   parameter meanings above are from the running text and the Part II dossier's description of the
   same scheme family.
3. Part II (`devleeschouwer2009_extraction.md`) cites this paper as "in press"; do not double-count
   the glucose column if Part II re-prints it.

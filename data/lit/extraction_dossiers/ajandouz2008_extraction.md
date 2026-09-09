# Ajandouz et al. 2008 — EXTRACTION (glucose / BSA / casein / glucose-protein at 60-100 C, pH 8.0 and 9.7)

**Source on disk:** `data/articles/ajandouz2008.pdf` (9 pp., born-digital Elsevier PDF; text layer clean,
Tables 1-3 verified against a 150-dpi raster of p. 3). Figures 1, 2 and 4 hold the time courses; only
Fig. 1 (glucose loss, pH 8.0) was digitised here (±3 % of axis). Read-only extraction, 2026-09-07.

## 0. Identity

| field | value |
|---|---|
| Title | "Effects of temperature and pH on the kinetics of caramelisation, protein cross-linking and Maillard reactions in aqueous model systems" |
| Authors | El Hassan Ajandouz, Véronique Desseaux, Sanaa Tazi, Antoine Puigserver (Université Paul Cézanne-Aix-Marseille III) |
| Venue | Food Chemistry 107 (2008) 1244-1252; received 19 Apr 2007, revised 15 Jun 2007, accepted 24 Sep 2007 |
| DOI as printed | `doi:10.1016/j.foodchem.2007.09.062` |

## 1. Why it matters to the model

A third laboratory measuring **glucose disappearance in water at 60-100 C**, with and without protein, and
the only paper in the corpus that runs the same charge at **two alkaline pH values (8.0 and 9.7)**. It
gives activation energies for glucose loss (~90 kJ/mol), for amino-group loss (92-143 kJ/mol) and for
UV-294 / A420 development (92-164 kJ/mol), and the within-study pH ratios (glucose loss 3.5-6x faster,
browning 10-50x faster at pH 9.7 than 8.0). It brackets the trunk's Martins 2005 glucose/glycine
constants (pH 6.8, 80-120 C) from the alkaline side: it tells the model how the glucose-side rates and
their Ea move when pH goes up, and confirms that glucose loss (caramelisation) has a lower Ea than
amine loss. Nothing here is a per-step constant; there is no reaction scheme.

## 2. Methods as they matter to a model

- **Charges:** glucose 0.2 M; protein (BSA or casein, Fluka) 5 mg/mL; or glucose 0.2 M + protein 5 mg/mL.
  Buffers: 0.2 M sodium phosphate for pH 8.0; 0.2 M sodium borate for pH 9.7. Final pH measured after
  heating (values not printed).
- **Vessel:** screw-cap tubes (volume not stated); water bath at 60, 70, 80, 100 C; ice quench.
- **Times:** pH 8.0: 240 min at 60 C, 180 min at 70 C, 120 min at 80 C, 40 min at 100 C. pH 9.7: 100 min
  at 60/70/80 C, 10 min at 100 C. Fig. 1 sampling at pH 8.0 / 100 C: 10, 20, 30, 40 min.
- **Replicates:** "Each kinetic point was measured at least in duplicate ... The mean deviation in the
  experimental data did not exceed 12%."
- **Analytics:** glucose by HPAEC-PAD (Dionex CarboPac PA-100, 5 mM NaOAc / 0.1 M NaOH, 0-250 pmol
  calibration); free amino groups by TNBS (Fields 1971) at 420 nm; A294 and A420 on diluted samples
  (< 1.5 AU); SDS-PAGE (10 % / 12 %) for cross-linking.
- **Rate constants:** "deduced from the slopes of the curves in the initial stages (r2 > 0.9)" — i.e.
  initial (zero-order-like) slopes, in % of initial per min; **the rate constants themselves are not
  tabulated anywhere**, only the Arrhenius Ea from them ("correlation coefficients higher than 0.98").
- **Quantification basis:** % of initial glucose / amino groups; absorbance units.

## 3. Tables (verbatim) and figure data

**Table 1.** "Activation energies associated with the disappearance of glucose, in heated glucose and
glucose-protein aqueous model systems at pH 8.0 and pH 9.7" — Activation energy (kJ per mol):

| | Glucose | Glucose + casein | Glucose + BSA |
|---|---:|---:|---:|
| pH 8.0 | 93 | 85 | 87 |
| pH 9.7 | 96 | 93 | 92 |

**Table 2.** "Activation energy Ea associated with the loss of free amino groups in protein-reducing sugar
model systems and in some kinds of food" (present-study rows; aw 1.0, 60-100 C, TNBS):

| System | Ea (kJ per mol) |
|---|---:|
| Casein pH 8.0 (alone) | 143 |
| BSA pH 8.0 (alone) | 116 |
| Glucose + casein pH 8.0 | 106 |
| Glucose + BSA pH 8.0 | 102 |
| Glucose + casein pH 9.7 | 105 |
| Glucose + BSA pH 9.7 | 92 |

Literature rows in Table 2 (not this study's data): casein pH 6.7, aw 0.7, 0-70 C, Van Slyke, 121 (Lea &
Hannan 1949); soy proteins aw 0.33-0.93, 30-80 C, 119 (Jokinen 1976); soy 6-18 % H2O, 80-130 C, 147
(Thompson 1976); ribose + BSA 14.1 % H2O, 85-145 C, 130 (Carpenter 1962); lactose + casein aw 0.33-0.98,
37-60 C, OPA, 156-117 (Malec 2002); herring meal 113; NFDM 46 and 39; pasta 31; rice + lysine 52.

**Table 3.** "Activation energy (Ea) associated with the browning process in glucose, glucose-protein and
milk-based systems^a" (A420, 60-100 C for present study; a: "All the milk-based systems were formed in
aqueous solution, pH 6.7"):

| System | Ea (kJ mol^-1) |
|---|---:|
| Glucose pH 8.0 | 164 |
| Glucose pH 9.7 | 126 |
| Glucose + casein pH 8.0 | 120 |
| Glucose + casein pH 9.7 | 92 |
| Glucose + BSA pH 8.0 | 130 |
| Glucose + BSA pH 9.7 | 95 |
| Lactose + casein, 110-150 C (Morales & Van Boekel 1998) | 125 |
| Glucose + casein, 110-150 C (Brands & Van Boekel 2001) | 121 |
| Lactose + casein, 90-130 C, kinetic modelling (Brands & Van Boekel 2002) | 71-159 |

**UV (A294) Ea, text only:** "152, 128 and 129 kJ per mol at pH 8.0, and 123, 107 and 101 kJ per mol at
pH 9.7, in the case of glucose, glucose-casein and glucose-BSA, respectively."

**Fig. 1** — "Kinetics of glucose loss in heated aqueous glucose (Glc) and glucose-casein model systems at
pH 8.0, at temperatures ranging from 60 C to 100 C." Digitised (% glucose remaining):

| T (C) | glucose alone | glucose + casein |
|---|---|---|
| 60 | 99 (30 min), 98 (60), 97 (120), 96 (180), 92 (240) | 96 (30), 97 (60), 94 (120), 90 (240) |
| 70 | 91 (30), 90 (60), 87 (120), 82 (180) | 96 (30), 89 (60), 86 (120), 82 (180) |
| 80 | 87 (30), 78 (60), 65 (90), 65 (120) | 92 (30), 83 (60), 77 (90), 68 (120) |
| 100 | 79 (10), 75 (20), 58 (30), 54 (40) | 78 (10), 70 (20), 57 (30), 55 (40) |

Derived initial first-order rates [D], glucose alone, pH 8.0 (from the digitised points, ±15 %):
0.00035 /min (60 C), 0.0011 (70 C), 0.0048 (80 C), 0.0154 (100 C); Arrhenius slope through these four gives
Ea ≈ 98 kJ/mol — consistent with the printed 93, so the digitisation is coherent. Zero-order equivalents
at 0.2 M: 0.07, 0.2, 0.8, 2.3 mmol/L/min.

Text ratios: at pH 9.7 glucose disappeared "at rates which were 3.5- to 6-fold higher than at pH 8.0, and
the temperature was found to have very little effect"; glucose-alone rates "equal to or higher than" with
protein; amino-group loss "2- to 3-fold and 2- to 16-fold higher in BSA and casein, respectively" with
glucose than without, "the lower the temperature, the higher the enhancing effects"; BSA alone "4-fold
more reactive than casein at 60 C" but equal at 100 C; UV and browning rate constants at pH 9.7 "10- to
50-fold higher than those measured at pH 8.0"; "UV absorbance of glucose accounts for 25-80 % of that of
glucose-protein mixtures, although the browning accounts for only 7-55 %", both rising with temperature;
Fig. 4 axis maxima at pH 8.0: A294 up to ~12 (glucose + BSA, 100 C, 40 min) and A420 up to ~1.2.

No reaction scheme was fitted; no multiresponse model.

## 4. What the repo could take

FIT-eligible (measured Ea, measured rates at four temperatures — but rates only via digitisation):

| candidate | numbers | note |
|---|---|---|
| Ea, glucose loss, water, no amine, pH 8.0 phosphate, 60-100 C | 93 kJ/mol (pH 9.7 borate: 96) | trunk k_glc_fru Ea 122.6 (Martins, pH 6.8, with glycine) vs k10 236 ± 63; Ajandouz's whole-glucose loss Ea of ~90 is a cross-lab check on the sugar-only limb |
| Ea, glucose loss with casein / BSA, pH 8.0 | 85 / 87 kJ/mol | amine lowers Ea by < 10 % |
| k_loss(glucose), pH 8.0, 100 C [D] | 0.015 /min (zero-order 2.3 mM/min) | trunk at 100 C, pH 6.8, 0.2 M glycine: k_glc_fru 1.6e-3 + k_schiff x 200 mM = 3.2e-3 + k10 4.4e-5 ≈ 0.005 /min total glucose loss — Ajandouz at pH 8.0 without any amine is 3x faster: pH 6.8 -> 8.0 roughly triples sugar-only glucose loss |
| k_loss(glucose), pH 8.0, 60/70/80 C [D] | 0.00035 / 0.0011 / 0.0048 /min | digitised; ±15 % |
| Ea, browning A420, glucose alone pH 8.0 / 9.7 | 164 / 126 kJ/mol | vs trunk Mel step 95.2 (Martins) — caramelisation browning is steeper in T than Maillard browning |
| Ea, browning, glucose + casein pH 8.0 / 9.7 | 120 / 92 kJ/mol | 120 matches Morales 1998 / Brands lactose-casein at pH 6.7 ("pH has little effect in the 6.7-8.0 range") |
| Ea, amino-group loss, glucose + casein / BSA pH 8.0 | 106 / 102 kJ/mol | vs trunk k_schiff Ea 96.8 |

Within-study ratios: pH 9.7 / pH 8.0 glucose loss 3.5-6x; pH 9.7 / pH 8.0 browning and A294 rates 10-50x;
Ea(pH 9.7)/Ea(pH 8.0) for browning 0.73-0.77 ("20-30 % decrease"); glucose-alone / glucose-protein browning
Ea 1.2-1.3 ("one fifth to one fourth higher"); amino loss with / without glucose 2-3x (BSA), 2-16x (casein);
caramelisation share of UV 25-80 %, of A420 7-55 % (rising with T).

Directional claims: glucose loss is unaffected or slowed by protein; amine loss is not increased from pH 8.0
to 9.7 (lysine still protonated) while glucose loss and browning are; Ea of glucose loss (~90) < Ea of amine
loss (~100) < Ea of browning (120-164 at pH 8.0) — "caramelisation occurs readily during the initial
stages"; protein cross-linking (SDS-PAGE) is faster with glucose at pH 8.0 but suppressed at pH 9.7.

## 5. Caveats

- **No rate constants printed**; only Ea. Rates are "initial slopes" in %/min on curves that are
  "variably linear", so they are neither first- nor zero-order constants; the digitised k above are my
  first-order approximations from Fig. 1 (pH 8.0 only; the pH 9.7 curves are "not shown").
- Buffers differ between the two pH levels (phosphate vs borate; borate complexes sugars), so the pH ratios
  carry a buffer effect.
- Proteins (BSA, casein) at 5 mg/mL — amino-group concentration ~0.3-0.4 mM, i.e. glucose is in ~500x
  excess; not a free amino acid; TNBS also counts arginine.
- Temperature range 60-100 C, below the trunk's 100-120 C cook window; 100 C runs are short (40 / 10 min).
- Final pH not reported; mean deviation up to 12 %.
- Table 3 attributes Ea 121 (110-150 C) to Brands & Van Boekel 2001, which contains no Ea (see
  `brands2001_extraction.md`); the value belongs to another Wageningen paper.
- Browning Ea are on raw A420 rates, not on a melanoidin concentration; A294 tracks "intermediate"
  products without identification.

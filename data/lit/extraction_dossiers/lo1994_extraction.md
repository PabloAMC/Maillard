# Lo, Westwood, McLellan, Selwood & Thornalley 1994 — EXTRACTION
### *"Binding and Modification of Proteins by Methylglyoxal under Physiological Conditions. A Kinetic and Mechanistic Study with Nα-Acetylarginine, Nα-Acetylcysteine, and Nα-Acetyllysine, and Bovine Serum Albumin"*
### Methylglyoxal (MGO) + Nα-acetyl-Arg / -Cys / -Lys (80 mM each, 1H NMR; stopped-flow UV for Cys) and + BSA ([14C]MGO), 100 mM Na-phosphate, pH 7.0-7.4, 37 °C. Rate and equilibrium constants for MGO-amine glycosylamines, the Arg dihydroxyimidazolidine / imidazolone, and the Cys hemithioacetal. **No pyrazines, no Strecker degradation (every alpha-amine is acetylated), no labels on ring carbons.** Bears on the pyrazine question only through the MGO sinks that compete with the Strecker step.

**Read 2026-10-31 from `data/articles/lo1994.pdf`** (7 pp., J. Biol. Chem., open access; publisher scan with a poor
OCR layer that drops every exponent). **All rate constants were read from 130 / 300 dpi renders of pp. 32301-32303**,
not from the text layer. Figs. 2 and 4-7 not digitised.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Theodore W. C. Lo, Marie E. Westwood, Antony C. McLellan, Trevor Selwood, Paul J. Thornalley (University of Essex) |
| Venue | *J. Biol. Chem.* **269** (51), 32299-32305, 23 Dec 1994; received 28 Jul 1994 |
| DOI | not printed. Crossref 2026-10-31: **10.1016/S0021-9258(18)31635-1** (same title, vol. 269, pp. 32299-32305) |

## 1. What this paper answers for the model

| need | what Lo 1994 gives | tag |
|---|---|---|
| pyrazine ratios for R | none | — |
| 13C / 15N label shares | none ([14C]MGO only to count protein-bound MGO) | — |
| aminoketone condensation k | none; no aminoketone can form (Nα-acetyl) | — |
| pyrrole vs pyrazine | none | — |
| peptide vs free amino acid | indirect: relative reactivity toward MGO (vs BSA = 1): aminoguanidine 0.178, L-arginine 0.0044, **L-lysine 0.0015** (pH 7.0, 37 °C) | [M] |
| **MGO sink: cysteine** | Nα-acetyl-Cys + MGO → hemithioacetal, >99 % conversion at 80 mM; **k2 = 4.1 ± 0.9 x 10⁴ M⁻¹ s⁻¹, k-2 = 7.5 ± 3.3 x 10⁻³ s⁻¹, K2 = 5.4 x 10⁶ M⁻¹** (unhydrated MGO); overall K (all MGO forms) **7.3 x 10⁴ M⁻¹**; equilibrium within ~2 s | [M] |
| **MGO sink: lysine epsilon-amine** | glycosylamine and bis-glycosylamine, reversible: **k1 6.8 x 10⁻³ M⁻¹ s⁻¹, k-1 5.0 x 10⁻⁴ s⁻¹, K1 13.6 M⁻¹; k2 1.9 x 10⁻² M⁻¹ s⁻¹, k-2 2.0 x 10⁻³ s⁻¹, K2 9.5 M⁻¹** (Model 3); then brown fluorescent oligomers (not characterised) | [M] |
| MGO sink: arginine | glycosylamine ⇌ dihydroxyimidazolidine: **k1 8.5 x 10⁻³ M⁻¹ s⁻¹, k-1 4.4 x 10⁻⁶ s⁻¹, k2 6.7 x 10⁻⁴ s⁻¹, k-2 "1.9⁻⁴" s⁻¹ (printed so; 1.9 x 10⁻⁴ by K2) [!], K1 1.9 x 10³ M⁻¹, K2 3.5**; imidazolone formation **2.4 x 10⁻⁶ M⁻¹ s⁻¹** (pH 7.0, 37 °C, n 40) | [M] |
| BSA | after 4 days, 55 % of MGO irreversibly bound, 30 % reversibly, 14 % degradation products, 1 % free; initial-rate law orders 0.88 (MGO) and 0.79 (BSA), k_obs 7.2 ± 0.3 x 10⁻³ "M⁻¹ s⁻¹" (unit inconsistent with fractional orders [!]) | [M] |
| Schiff base | not detected in solution (acetal 1H ~8.5 ppm appears only after lyophilisation) [Q] | [Q] |

**Answers.**

1. **No pyrazine number; it does not touch P1x / P1y vs P2b.**
2. **Cysteine sequesters MGO fast and almost completely at 37 °C** [M]: with K ≈ 7.3 x 10⁴ M⁻¹, 20 mM free thiol leaves
   ~0.07 % of MGO unbound at equilibrium [Z]. If anything like this holds at 100-130 °C (K will fall with
   temperature; not measured here), the MGO available to the Strecker step in a cysteine pot is far lower than total
   MGO. That is a candidate cause for every Ho-network cysteine row sitting below R = 4 and for lower DMP in Cys pots
   (EX-052 red team sec. 5) — a hypothesis, not a test.
3. **Reversible MGO-amine adducts dominate over Schiff bases** at pH 7, 37 °C [Q]; the epsilon-amine binds MGO weakly
   (K1 13.6 M⁻¹) — a reservoir, not a sink, at the concentrations of model pots.

## 2. Conditions [M]

| item | value |
|---|---|
| Reactants | MGO (from 1,1-dimethoxypropanone, distilled; assayed enzymatically) 80 mM + Nα-acetyl amino acid 80 mM (NMR); 0.25-2.0 mM equimolar (stopped flow, Cys); [14C]MGO 1.0 mM + BSA 0.10 mM (protein); scavenger competition 0-29 mM |
| Buffer / pH | Na-phosphate 100 mM (200 mM for the imidazolone preparation); **pH 7.0 in the Methods and figure legends, pH 7.4 in the text and abstract** [!] |
| Temperature | 37 °C; NMR samples in 50 % D2O |
| O2 | air, sterile; the imidazolone forms by autoxidation of a 1,5-dihydroimidazolone [P] |
| Analysis | 1H NMR (270 MHz) time courses, fitted with the GIT Gear integrator; stopped-flow UV 288 nm (hemithioacetal ε 98 ± 5 M⁻¹ cm⁻¹); RP-HPLC 320 nm (imidazolone, 2-methylquinoxaline IS); ultrafiltration / perchloric precipitation for BSA |

## 3. Quantitative results

All in sec. 1 table. Model 2 assumed (from ref 28) k1 3.1 x 10⁻³ s⁻¹, k-1 0.16 s⁻¹, K1 0.019 for MGO hydrate
dehydration [C]. Fig. 2's y-axis is labelled "mM" with values 0-0.10 for an 80 mM pot [!] (axis in M). HSA estimates
(46 mg/mL, 80 nM MGO): reversible binding Arg : Lys : Cys ≈ 110 : 51 : 0.6 relative to free MGO; >99 % of plasma MGO
protein-bound [M]. (The authors' in vivo modification-rate estimate on p. 32304 is garbled in the OCR and was not
re-read from a render; not used.)

## 4. Mechanism evidence and its kind

- **Measured intermediates** (kind 4), Essex (Thornalley): glycosylamines, dihydroxyimidazolidine, hemithioacetal
  identified by 1H NMR; the imidazolone isolated and characterised (1H / 13C NMR, FAB-MS m/z 269 M+1).
- **Fit** (kind 7) for the rate constants (Models 1-3, simplified schemes).
- **Competition kinetics** with scavengers (kind 3, trapping).

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 / CF-077 / P1x-P1y-P2b | LEAVES UNCHANGED | no pyrazine data |
| the Ho-network cysteine rows (R < 4, low DMP) | offers a mechanism [Z] | fast reversible MGO-thiol binding cuts free MGO |
| CF-073 (Strecker yield from MGO, Weenen vs Hofmann) | LEAVES UNCHANGED | |
| 3-DX / dicarbonyl trapping by cysteine (evidence map 2.2) | SUPPORTS for MGO (37 °C) | a measured carbonyl-thiol K and rate; MGO, not 3-DX |
| protein-bound amine (triage 9) | adds MGO-lysine / arginine constants at 37 °C | far from process temperature |

## 6. Disagreements found

1. **pH 7.0 vs 7.4 within the paper** (Methods / legends vs text / abstract). Cause: unknown; constants are
   probably pH 7.0 (the experimental descriptions).
2. **k-2 printed "1.9⁻⁴ s⁻¹"** (Model 1); K2 = 3.5 = 6.7 x 10⁻⁴ / 1.9 x 10⁻⁴ confirms 1.9 x 10⁻⁴.
3. No conflict with other dossiers found; the MGO-cysteine K is not keyed anywhere yet.

## 7. Ledger-ready measurements

1. MGO + Nα-acetyl-Cys: k2, k-2, K2, overall K 7.3 x 10⁴ M⁻¹ (pH 7.0, 37 °C) — the carbonyl-thiol equilibrium for the
   cysteine-pot MGO budget (temperature transfer not licensed).
2. MGO + Nα-acetyl-Lys: Model 3 constants — MGO reservoir on epsilon-amines.
3. Relative MGO reactivities (aminoguanidine 0.178, Arg 0.0044, Lys 0.0015 vs BSA).

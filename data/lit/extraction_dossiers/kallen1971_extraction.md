# Kallen 1971 — EXTRACTION (cysteine + formaldehyde → thiazolidine-4-carboxylate (TC) equilibria at 25 °C, ionic strength 1.0 M (KCl), pH 5.57 / 6.56 / 7.45 in 0.01 M phosphate; per-species formation constants for TC from CYS⁺ / CYS± / CYS⁻; N-hydroxymethyl and hemithioacetal adduct constants; the full microscopic pKa set of cysteine; one approach-to-equilibrium figure, NO rate constants)

### THE CLEANEST PER-SPECIES EQUILIBRIUM DECOMPOSITION OF A CYSTEINE THIAZOLIDINE IN THE CORPUS: the thiazolidine forms ~217× more favourably from the cysteine monoanion (CYS⁻) than from the zwitterion (CYS±), so the apparent K rises ~20× between pH 5.6 and 7.5. But the carbonyl is formaldehyde (hydrated, no ring-chain tautomerism), the temperature is 25 °C only, and the rates are in the companion paper (Kallen 1971, JACS 93, 6236; ref 9), which is NOT on disk.

**Source on disk:** `data/articles/kallen1971.pdf` (9 pp., JACS pp. 6227-6235; p. 6227 starts with the tail of an
unrelated silane paper). Prose and Tables I, IV re-typed from the `pdftotext -layout` text layer. **The text layer
drops the row labels of Table II and garbles Table III's exponents**; both were re-typed from 130-400 dpi page
renders (scratchpad `img/kal-5.png`, `kal-6.png`, `kal6crop-6.png`). Figure 4 read from a 300 dpi crop
(`img/kalfig4-6.png`); every number read off it is tagged **[D]**.

## 0. Identity

| field | value |
|---|---|
| Title | "Equilibria for the Reaction of Cysteine and Derivatives with Formaldehyde and Protons" |
| Authors | Roland G. Kallen — Department of Biochemistry, School of Medicine, University of Pennsylvania |
| Venue | Journal of the American Chemical Society **93**(23), 6227-6235 (November 17, 1971); received August 3, 1970 |
| DOI as printed | **not printed** |
| Companion | ref 9: R. G. Kallen, J. Amer. Chem. Soc. 93, 6236 (1971) — the mechanism / pH-rate paper for TC formation ("exhibits a complex pH-rate profile"). **Not on disk; this is the paper that would carry forward/reverse rate constants.** |
| Paper type | Physical-organic equilibrium study (titrimetry, UV spectrophotometric titration, DTNB thiol assay, formol titration) |

## 1. Why it matters (B59b: Cys + pentose ⇌ TTCA reservoir)

The engine needs K, k_f, k_r for Cys + PENT ⇌ TTCA and their pH and T dependence. The corpus has Cao 2017
(`cao2017_extraction.md`, Cys + Xyl, 90 °C, K ≈ 60 L/mol). Kallen delivers, for the formaldehyde analogue:

1. **Per-microspecies formation constants** (Table III): which cysteine species the thiazolidine equilibrium couples to,
   and how strongly. This is the transferable content: the **pH shape** of K_app (eq 5), not its magnitude.
2. **Apparent K vs pH** at 5.57, 6.56, 7.45 (Table IV): 6.53 × 10⁵ → 2.32 × 10⁶ → 1.30 × 10⁷ M⁻¹.
3. **The TC ring-nitrogen pKa = 6.24** (Table I) — for TTCA the analogous pKa sets where the reservoir's K_app turns over.
4. **No rate constants.** Only Fig. 4 (approach to equilibrium from TC hydrolysis, two pH values) and a remark that the
   observed relaxation rates agree with ref 9.

**Delivers for B59b:** the pH-shape of K and the reacting species — **yes**; a constant transferable to an aldose — **no**
(formaldehyde ≠ aldose; 25 °C only; I = 1.0 M); rates — **no** (only a crude relaxation time, arithmetic 3).

## 2. Methods as they matter to a model

- **Solvent / ionic strength:** deionised water, KCl to I = 1.0 M throughout. 25 °C (± 0.1 °C for spectrophotometry);
  ΔH for carbinolamines from 10 / 25 / 40 °C.
- **TC equilibrium (Table IV):** 50 mL reaction mixture, **TC 0.5-7.5 mmol/L** initial (synthesised TC, mp 189-191 °C),
  0.001 M EDTA (Methods) / "0.01 M EDTA" (Table IV and Fig. 4 footnotes — **inconsistent in the source**), 0.01 M
  potassium phosphate at **pH 5.56(5.57), 6.56, 7.45**, 25 °C, protected from light. Aliquots (2.7 mL) quenched into
  0.3 mL 1.0 M phosphate (90 % free base) + 0.02 mL 0.01 M DTNB; free thiol by ΔA412 (ε 13,600 M⁻¹ cm⁻¹) against a
  cysteine standard. Equilibrium approached from TC hydrolysis; at pH 5.73 also from the synthesis side
  (K_app ≈ 7 × 10⁵ M⁻¹, average). Free formaldehyde hydrate at equilibrium = free cysteine (1:1 stoichiometry).
- **Carbinolamine constants:** formol titration of 0.01 M S-methylcysteine (SMC) buffers (10 % free base) vs free
  formaldehyde hydrate, 10 / 25 / 40 °C (Fig. 3, eq 4).
- **Hemithioacetal (K_H′, NAC):** spectrophotometric, from a manuscript submitted elsewhere (Kallen & Frederick).
- **Carbinolamine of cysteine itself (26.8 M⁻¹):** kinetic, pH 10-11, from ref 9.
- **Formaldehyde species:** all constants are against **formaldehyde hydrate** F (CH₂(OH)₂), water activity 1.0.
- Concentrations in mmol/L: TC 0.5-7.5 mmol/L; cysteine for titrations 10 mmol/L; spectrophotometric 0.10-0.15 mmol/L.

## 3. Tables re-typed

### Table I. Macroscopic proton dissociation constants at 25 °C (this work: I = 1.0 M; n in parentheses)

| compound | pG1′ | pG2′ | pG3′ | ref |
|---|---|---|---|---|
| S-methyl-L-cysteine | 2.22 (19) | 8.88 (28) | | this work |
| S-methyl-L-cysteine | n.d. | 8.75 | | 10e |
| S-ethyl-L-cysteine | 2.03 | 8.60 | | 3d (I 0.1-0.2 M) |
| N-acetyl-L-cysteine | 3.08 (29) | 9.51 (52) | | this work |
| L-cysteine methyl ester | 6.77 (51) | 8.99 (51) | | this work (Noyes method) |
| L-cysteine ethyl ester | 6.69 | 9.17 | | 10d (I ≈ 0.05 M) |
| L-cysteine betaine | n.d. | 8.65 | | 10e (I 0.15 M) |
| N-acetyl-S-benzyl-L-cysteine | 3.07 (17) | | | this work |
| **L-cysteine** | **2.00 (13)** | **8.38 (71)** | **10.37 (71)** | **this work** (Noyes; ranges pG2′ 8.34-8.44, pG3′ 10.30-10.45) |
| L-cysteine | n.d. | 8.30 | 10.40 | 10e |
| L-cysteine | 1.71 | 8.33 | 10.78 | 10b (I = 0) |
| L-cysteine | n.d. | 8.27 | 10.42 | 10g (I 0.1 M) |
| L-cysteine | n.d. | 8.07 | 9.95 | 10r (I 1.0 M) |
| **thiazolidine-4-carboxylate** | **1.63 (25)** | **6.24 (54)** | | **this work** |
| thiazolidine-4-carboxylate | 1.51 | 6.21 | | 3d |

### Table II. Microscopic pK′ of cysteine at 25 °C, I = 1.0 M (this work, columns d and e; selected comparison column f = ref 10p)

Notation (Scheme I): subscript 1 = COOH, 2 = SH, 3 = NH3⁺; pK12′ = thiol loss from the zwitterion (NH3⁺ retained);
pK13′ = ammonium loss with SH retained; pK123′ / pK132′ = the second group from those monoanions.

| constant | d (Wegscheider, SMC) | e (spectrophotometric) | f (10p) |
|---|---|---|---|
| pK1′ | 2.00 | 2.00 | 1.70 |
| pK2′ | 7.44 | 7.44 | 7.45 |
| pK3′ | 6.88 | 6.88 | 6.77 |
| **pK12′** (zwitterion → ⁻S, NH3⁺) | **8.55** | **8.54** | 8.53 |
| **pK13′** (zwitterion → SH, NH2) | **8.88** | **8.90** | 8.86 |
| pK21′ | 3.10 | 3.10 | 2.79 |
| pK23′ | 8.32 | 8.32 | 8.41 |
| pK31′ | 4.00 | 4.02 | 3.80 |
| pK32′ | 8.88 | 8.88 | 9.09 |
| **pK123′** (⁻S,NH3⁺ → ⁻S,NH2) | **10.20** | **10.22** | 10.36 |
| **pK132′** (SH,NH2 → ⁻S,NH2) | **9.87** | **9.85** | 10.03 |
| pK231′ | 4.99 | 4.99 | 4.74 |
| pK3-12′ | 3.67 | 3.66 | 3.47 |
| pK2-1′ | −4.88 | −4.88 | −5.06 |
| pK2-3′ | −0.56 | −0.56 | −0.68 |
| **pK12-13′** (tautomer ratio, log [⁻S,NH3⁺]/[SH,NH2]) | **0.33** | **0.36** | 0.33 |
| pK23-12′ | −5.22 | −5.23 | −5.62 |

Other literature columns (g-l) give pK12′ 8.45-8.66, pK13′ 8.31-8.86, pK123′ 9.58-10.45, pK132′ 9.73-10.51,
tautomer log ratio −0.15 to 0.35. NAC (this work): pK1′ 3.08, pK2′ 8.48 (estimated), pK12′ 9.51, pK21′ 4.11. CME: pK2′ 7.44,
pK3′ 6.88, pK23′ 8.32, pK32′ 8.88. Text: tautomerisation constant [⁻S,NH3⁺]/[SH,NH2] "is about 2".

### Table III. Equilibrium constants, I = 1.0 M, water activity 1.0 (F = formaldehyde hydrate)

| equilibrium | T, °C | SMC | NAC | CYS |
|---|---|---|---|---|
| L1′ = [RNHCH2OH][HOH]/([RNH2][F]) | 10 | 21.0 | | |
| | 25 | 15.9 | | 26.8 (kinetic, pH 10-11, ref 9) |
| | 40 | 11.8 | | |
| L3′ = [RN(CH2OH)2][HOH]/([RNHCH2OH][F]) | 10 | 1.54 | | |
| | 25 | 1.02 | | |
| | 40 | 0.69 | | |
| K_H′ = [RSCH2OH][HOH]/([RSH][F]) (hemithioacetal) | 25 | | 700 | |
| **K_F1′ = [TC⁺]/([CYS⁺][F])** | 25 | | | **2.24 × 10⁵** (from K_F2′·G1′/K_TC1′; 0.01 M phosphate, 0.01 M EDTA) |
| **K_F2′ = [TC±]/([CYS±][F])** | 25 | | | **5.25 ± 0.5 × 10⁵** (eq 5 fit) |
| **K_F3′ = [TC⁻]/([CYS⁻][F])** | 25 | | | **1.14 ± 0.01 × 10⁸** (eq 5 fit) |
| K_TF1′ = [HMTC⁻]/([TC⁻][F]) | 25 | | | 7.0 |
| K_TF2′ = [HMTC±]/([TC±][F]) | 25 | | | 0.7 |

Units M⁻¹ throughout. SMC pG2′ = 9.23, 8.88, 8.19 at 10, 25, 40 °C (footnote b). ΔH1 = −3.3 and ΔH3 = −4.8 kcal/mol
(van 't Hoff of L1′, L3′). **The "×10⁵" exponents on K_F1′ and K_F2′ are ambiguous in print (5 vs 6); they are fixed at
10⁵ by three internal checks** (arithmetic 1). CYS⁻ here is the macroscopic monoanion (both tautomers ⁻S,NH3⁺ and SH,NH2).

Eq 5 as printed (pH > 4): 1/K_Fapp′ = (1 − α_TC⁻)/K_F2′ + α_TC⁻ (1 + G3′/[H⁺])/K_F3′, α_TC⁻ = fraction of TC as
anion. **I could not make this printed form reproduce Table IV**; the Scheme III species bookkeeping
K_app = K_F2′·x(CYS±)/(1 − α_TC⁻) = K_F3′·x(CYS⁻)/α_TC⁻ does (arithmetic 1). Footnote 38: equivalently, for pH 5.5-7.5, a plot of K_Fapp′/α_CYS⁻ against α_CYS±/α_CYS⁻ has slope
K_F2′ and intercept K_F3′.

### Table IV. Apparent K for TC formation, 25 °C, I = 1.0 M, 0.01 M phosphate, 0.01 M EDTA

| initial TC (mmol/L) | Cys eq (10⁻⁵ M) pH 5.57 | 6.56 | 7.45 | TC eq (mmol/L) 5.57 | 6.56 | 7.45 | K_Fapp′ (10⁵ M⁻¹) 5.57 | 6.56 | 7.45 |
|---|---|---|---|---|---|---|---|---|---|
| 7.5 | 10.8 | 5.96 | 2.62 | 7.39 | 7.44 | 7.47 | 6.33 | 20.9 | 146 |
| 6.0 | 9.96 | 5.18 | 2.27 | 5.90 | 5.95 | 5.98 | 5.94 | 22.2 | 116 |
| 4.0 | 7.61 | 4.28 | 1.83 | 3.92 | 3.96 | 3.98 | 6.76 | 21.6 | 119 |
| 2.0 | 5.31 | 2.91 | 1.28 | 1.95 | 1.97 | 1.99 | 6.91 | 23.3 | 121 |
| 1.0 | 3.77 | 1.93 | 0.83 | 0.96 | 0.98 | 0.99 | 6.75 | 26.3 | 144 |
| 0.5 | 2.69 | 0.98 | 0.60 | 0.47 | 0.49 | 0.49 | 6.49 | 24.6 | 136 |
| **mean** | | | | | | | **6.53** | **23.2** | **130** |
| s.d. | | | | | | | 0.32 | 1.84 | 12.0 |

K_Fapp′ = [TC]eq / ([Cys]eq [F]eq) with [F]eq = [Cys]eq; constancy over a 15-fold TC range confirms 1:1 stoichiometry.

### Figure 4 (approach to equilibrium, TC hydrolysis, 25 °C, pH 5.39 and 6.19) [D]

ΔA412 at t = 0 ≈ 0.27-0.30; plateau ≈ 1.12 (pH 6.19) and ≈ 1.43 (pH 5.39). Half of the rise reached at
**≈ 0.4 h (pH 6.19)** and **≈ 1.2-1.4 h (pH 5.39)** [D, ±30 %]; plateau by ~2 h and ~10 h respectively.
Caption initial TC "7.5 × 10⁻³ M" (exponent faint), right axis "free thiol × 10⁵ (M)" plateau ≈ 0.8-1.0 × 10⁻⁵ M —
**about 10× below Table IV's free Cys at 7.5 mmol/L TC**; the source does not reconcile them (dilution in the assay
cuvette is only 2.7/3.02). Treat Fig. 4 as a timescale only.

## 4. Arithmetic (mine)

**1. Fixing the exponents and the eq 5 form.** Cycle: K_F3′/K_F2′ = K_TC2′/G2′ = 10^(8.38 − 6.24) = 138.
1.14 × 10⁸ / 138 = 8.3 × 10⁵ ≈ K_F2′ = 5.25 × 10⁵ (within 1.6×; would be 22× off with 10⁶). Text: "formation of TC± from
CYS± is more than two orders of magnitude less favorable than TC⁻ from CYS⁻": 1.14 × 10⁸ / 5.25 × 10⁵ = **217** ✔
(10⁶ would give 22). K_F1′ = K_F2′ · G1′/K_TC1′ = 5.25 × 10⁵ × 10^(−2.00 + 1.63) = **2.24 × 10⁵** ✔.
Species check of Table IV: K_app = K_F2′·x(CYS±)/(1 − α_TC⁻) at pH 5.57 [x(CYS±) ≈ 0.998, α_TC⁻ = 0.176] =
5.25 × 10⁵/0.824 = **6.4 × 10⁵** vs 6.53 × 10⁵ ✔; K_app = K_F3′·x(CYS⁻)/α_TC⁻ at pH 7.45 [x(CYS⁻) = 0.105, α_TC⁻ = 0.942]
= **1.27 × 10⁷** vs 1.30 × 10⁷ ✔.

**2. pH sensitivity.** log K_app rises 0.55 per unit pH between 5.57 and 6.56 and 0.83 per unit between 6.56 and 7.45;
×20 over 1.9 pH units. Mechanistic reading: K_app ∝ fraction of cysteine as monoanion / fraction of TC as anion; below
the TC ring-N pKa (6.24) the protonated ring is strongly destabilised.

**3. Relaxation rate (crude).** With [Cys]eq/[TC] ~ 1 %, the relaxation constant ≈ the hydrolysis (ring-opening +
release) constant: k_r ≈ ln 2 / t½ ≈ **5 × 10⁻⁴ s⁻¹ (pH 6.19)** and **≈ 1.5 × 10⁻⁴ s⁻¹ (pH 5.39)** at 25 °C [D, factor ~2].
Then k_f = K_app · k_r ≈ 1.5 × 10⁶ × 5 × 10⁻⁴ ≈ **~10³ M⁻¹ s⁻¹ (pH 6.2, formaldehyde, 25 °C)** — a Z-grade order of magnitude
only; the exact pH-rate law is in ref 9.

**4. Against Cao 2017.** Cys + xylose K ≈ 60 L/mol at 90 °C vs formaldehyde 6.5 × 10⁵ M⁻¹ (pH 5.57, 25 °C): **~10⁴× lower**.
Plausible contributors, none quantified here: aldose open-chain fraction (≪ 1 %), steric/electronic substitution at C2,
and temperature (ΔH of TC formation not measured; carbinolamine ΔH −3.3 kcal/mol would lower K by only ~2.5× from 25 to
90 °C). **Not a derivation — a consistency check that the magnitudes are not contradictory.**

## 5. What the repo could take

| item | value | conditions | grade |
|---|---|---|---|
| K_F3′ (TC⁻ from CYS⁻ + F) | 1.14 × 10⁸ M⁻¹ | 25 °C, I 1.0 M | **F** (eq 5 nonlinear fit to Table IV) |
| K_F2′ (TC± from CYS± + F) | 5.25 × 10⁵ M⁻¹ | same | **F** |
| K_F1′ (TC⁺ from CYS⁺ + F) | 2.24 × 10⁵ M⁻¹ | same | **F** (derived by author from K_F2′) |
| K_Fapp′ vs pH | 6.53 × 10⁵ / 2.32 × 10⁶ / 1.30 × 10⁷ M⁻¹ | pH 5.57 / 6.56 / 7.45, 25 °C | **M** |
| TC ring-N pKa | 6.24 (TC COOH 1.63) | 25 °C, I 1.0 M | **M** |
| cysteine macroscopic pKa | 2.00 / 8.38 / 10.37 | 25 °C, I 1.0 M | **M** |
| cysteine microscopic pKa (thiol from zwitterion 8.55; amine 8.88; tautomer ratio ≈ 2) | Table II | 25 °C, I 1.0 M | **M/F** (spectrophotometric + Wegscheider) |
| carbinolamine L1′ (SMC) | 15.9 M⁻¹, ΔH −3.3 kcal/mol | 10-40 °C | **M** |
| hemithioacetal K_H′ (NAC) | 700 M⁻¹ | 25 °C | **M** (cited manuscript) |
| relaxation k ≈ 1.5-5 × 10⁻⁴ s⁻¹; k_f ~10³ M⁻¹ s⁻¹ | pH 5.4-6.2, 25 °C | **Z** (figure read, arithmetic 3) |
| **directional check:** TTCA-type K_app should RISE with pH across 5-7.5 (≈ ×20 for formaldehyde) and the reservoir couples to the cysteine monoanion | | **Z** (shape transfer is an assumption) |
| ! Rate law and k's are in ref 9 (JACS 93, 6236), not here — the paper to fetch for B59b rates | | **!** |
| ! Formaldehyde is a hydrate-dominated, non-cyclic aldehyde; magnitudes do not transfer to ribose/xylose | | **!** |
| ! I = 1.0 M KCl and 25 °C only; no ΔH for thiazolidine formation | | **!** |
| ! Eq 5 as printed ("G3′") does not reproduce Table IV; use Scheme III bookkeeping | | **!** |
| ! EDTA 0.001 M (Methods) vs 0.01 M (Table IV, Fig. 4 footnotes) | | **!** |
| ! Abstract gives N,N-dihydroxymethyl constant "16.0 M⁻¹" (should be M⁻²; text uses L2′ = 16.2) and swaps "neutral/anionic" for the HMTC constants relative to Table III (Table III: HMTC⁻ 7.0, HMTC± 0.7) | | **!** |

## Use in the model

- **B59b, pH dependence:** use Kallen's species structure, not its numbers: write K_app(pH) = K_anion·x(Cys⁻)/α(TTCA⁻) +
  (zwitterion term), with the TTCA ring-N pKa (6.24 for TC; TTCA's own pKa should be sourced) controlling the turnover.
  Anchor the magnitude on Cao 2017 and let Kallen supply the slope (~×20 per 1.9 pH units at 25 °C) as a declared
  assumption with a wide band.
- **Rates:** only a Z-grade timescale (hours at 25 °C, mM, formaldehyde). At 90-140 °C and with Cao's observed 0.5 h
  near-completion, the fast-equilibrium assumption is consistent; this paper neither proves nor bounds it for aldoses.
- **Fetch:** Kallen 1971, JACS 93, 6236 (mechanism and pH-rate profile of TC formation).

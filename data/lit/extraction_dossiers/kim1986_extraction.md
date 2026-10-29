# Kim, Yun & Park 1986 — EXTRACTION (cysteine + cinnamaldehyde, 25 °C, I = 0.1 M NaCl, pH 1-13: Michael addition of cysteine to the C=C (step 1) then thiazolidine formation of the monoadduct's aldehyde with a SECOND cysteine (step 2); pH-rate profile of step 1, a pH 8.0 time course, k1 / k−1 / k2 / k−2 at 13-38 °C, activation parameters)

### STEP 2 IS THE TRANSFERABLE PART: a saturated aliphatic aldehyde (the S-cysteinyl-3-phenylpropanal monoadduct) + cysteine → thiazolidine, k2 = 122 M⁻¹ s⁻¹, K2 = 1.49 × 10⁶ M⁻¹ at pH 8.0, 25 °C, Ea 6.1 kcal/mol — but at ONE pH only. The full pH-rate profile (Tables 1-2, eq 10) is for the Michael addition, which no aldose does.

**Source on disk:** `data/articles/Kim1986.pdf` (5 pp.: tail of the preceding paper on p. 25, the article pp. 25-28,
head of the next paper on p. 29). **The text layer is badly garbled (Korean OCR overlay; exponents lost)**; every table
and equation below was re-typed from 140 dpi page renders and 300 dpi crops (scratchpad `img/kim-2.png`, `kim-3.png`,
`kim-4.png`, `kimT3-3.png`, `kimT4-4.png`, `kimT4b-4.png`, `kimEq9-3.png`, `kimEq10-4.png`). Figures 1-3 not digitised.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics and Thermodynamic Studies on the Reaction of Cysteine with Cinnamaldehyde" (running head: "The Reaction mechanism of Cysteine with Cinnamaldehyde") |
| Authors | Tae-Rin Kim (Department of Chemistry, Korea University, Seoul); Se-Joon Yun, Byung Bin Park (Kongju National Teachers College) |
| Venue | Bulletin of the Korean Chemical Society **7**(1), 25-28 (1986); received July 8, 1985 |
| DOI as printed | **not printed** |
| Paper type | Kinetics by UV (cinnamaldehyde at 293 nm) and polarography (DME vs SCE: cysteine wave −0.25 V, cinnamaldehyde −1.40 V) |

## 1. Why it matters (B59b)

The engine needs the forward/reverse rates of Cys + carbonyl ⇌ thiazolidine and their pH, T dependence. This paper
gives:

1. **A thiazolidine-forming step with an aliphatic aldehyde, rates in both directions, at four temperatures** (Table 4,
   step 2). Only at pH 8.0.
2. **A pH-rate profile (pH 1-13) for the thiol Michael addition** (step 1) with a microscopic-species rate law (eqs 9-10).
   This identifies thiolate as the reactive nucleophile for **C=C addition**, which is not the TTCA reaction.
3. **A pH 8.0 time course** (Table 3) of cinnamaldehyde, cysteine, mono- and di-adduct — the only data the four rate
   constants were fitted to.

**Delivers for B59b:** an order-of-magnitude rate and K for an aliphatic aldehyde + cysteine → thiazolidine at pH 8, 25 °C,
with a small Ea — **partially**; a pH-rate law transferable to an aldose — **no** (the profile belongs to step 1).

## 2. Methods as they matter to a model

- **Temperature** 25 °C (rate constants also at 13, 32, 38 °C, pH 8.0). **Ionic strength 0.1** (NaCl). N₂ atmosphere.
- **Buffers (Table 1):** HCl + KCl (pH 1-3.5), HOAc + NaOAc (4.0-5.5), KH₂PO₄ + borax (6.0-9.0), NaOH + NaCl (9.5-12.0).
  Buffer concentrations not stated; "general base effect was actually negligible".
- **Spectrophotometric runs:** 98 mL buffer + 1 mL of 2 × 10⁻³ M cinnamaldehyde in **alcohol** + 1 mL cysteine solution
  → **cinnamaldehyde 0.020 mmol/L**, **cysteine 0.10-1.00 mmol/L** (pseudo-first order, 1 % ethanol). k1 from the slope
  of k_obs vs [Cys] (Fig. 1, pH 8.5).
- **Polarographic runs:** 5 mL 1 × 10⁻³ M cinnamaldehyde + 10 mL 1 × 10⁻³ M cysteine in buffer → **cinnamaldehyde
  0.33 mmol/L, cysteine 0.66 mmol/L**. Mono- and di-adduct inferred by mass balance (eqs 3-6), not measured directly.
- **Fit:** eq 7 solved "by simultaneous equation" on Fig. 3 data: −d[Cys]/dt = k1[CA][Cys] − k−1[CACys1] + k2[Cys][CACys1]
  − k−2[CA(Cys)2]. No uncertainties given.
- Cysteine pKa taken from Wrathall, Izatt & Christensen 1964 (ref 10): K1 = 4.07 × 10⁻⁹ (pK 8.39), K2 = 6.98 × 10⁻¹¹
  (pK 10.16). Kallen 1971 (`kallen1971_extraction.md`) says several of that reference's values disagree with other work.

## 3. Tables re-typed

### Table 1. k1 by spectrophotometry, 25 °C, I = 0.1 (title prints "Cystine" — typo for cysteine)

| pH | buffer | k1 obs (M⁻¹ s⁻¹) | log k1 + 3 | k1 calc (eq 10) | log k1 + 3 calc |
|---|---|---|---|---|---|
| 1.0 | HCl + KCl | — | — | 5.50 × 10⁻³ | 0.74 |
| 2.0 | HCl + KCl | 5.62 × 10⁻³ | 0.75 | 5.50 × 10⁻³ | 0.74 |
| 2.5 | HCl + KCl | 5.75 × 10⁻³ | 0.76 | 5.51 × 10⁻³ | 0.74 |
| 3.0 | HCl + KCl | 5.85 × 10⁻³ | 0.77 | 5.53 × 10⁻³ | 0.74 |
| 3.5 | HCl + KCl | 5.59 × 10⁻³ | 0.78 | 5.61 × 10⁻³ | 0.75 |
| 4.0 | HOAc + NaOAc | 6.00 × 10⁻³ | 0.78 | 5.84 × 10⁻³ | 0.77 |
| 4.5 | HOAc + NaOAc | 7.02 × 10⁻³ | 0.85 | 6.58 × 10⁻³ | 0.82 |
| 5.0 | HOAc + NaOAc | 1.02 × 10⁻² | 1.08 | 8.92 × 10⁻³ | 0.95 |
| 5.5 | HOAc + NaOAc | 1.62 × 10⁻² | 1.21 | 1.63 × 10⁻² | 1.21 |
| 6.0 | KH₂PO₄ + borax | 5.72 × 10⁻² | 1.76 | 3.95 × 10⁻² | 1.60 |
| 6.5 | KH₂PO₄ + borax | 1.39 × 10⁻¹ | 2.14 | 1.12 × 10⁻¹ | 2.05 |
| 7.0 | KH₂PO₄ + borax | 4.88 × 10⁻¹ | 2.69 | 3.34 × 10⁻¹ | 2.52 |
| 7.5 | KH₂PO₄ + borax | 1.00 | 3.00 | 9.62 × 10⁻¹ | 2.98 |
| 8.0 | KH₂PO₄ + borax | 3.52 | 3.55 | 2.43 | 3.39 |
| 8.5 | KH₂PO₄ + borax | 3.79 | 3.58 | 4.67 | 3.67 |
| 9.0 | KH₂PO₄ + borax | 4.98 | 3.70 | 6.39 | 3.81 |
| 9.5 | NaOH + NaCl | 4.75 | 3.68 | 6.48 | 3.81 |
| 10.0 | NaOH + NaCl | 4.40 | 3.64 | 4.90 | 3.69 |
| 10.5 | NaOH + NaCl | 1.25 | 3.10 | 2.65 | 3.42 |
| 11.0 | NaOH + NaCl | 0.90 | 2.77 | 1.10 | 3.04 |
| 11.5 | NaOH + NaCl | 0.29 | 2.47 | 4.12 × 10⁻¹ | 2.61 |
| 12.0 | NaOH + NaCl | — | — | 1.68 × 10⁻¹ | 2.23 |
| 12.5 | — | — | — | 8.77 × 10⁻² | 1.94 |
| 13.0 | — | — | — | 6.19 × 10⁻² | 1.79 |

(pH 13.5 and 14.0 rows empty in print. "log k1 + 3" for pH 3.5-4.0 observed rows do not match the k1 values exactly
(5.59 × 10⁻³ → 0.75, printed 0.78) — as printed.)

### Table 2. k1 by polarography, 25 °C, I = 0.1

| pH | 6.0 | 6.5 | 7.0 | 7.5 | 8.0 | 8.5 | 9.0 | 9.5 | 10.0 | 10.5 | 11.0 | 11.5 | 12.0 | 12.5 | 13.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k1 obs (M⁻¹ s⁻¹) | 3.96 × 10⁻² | 1.02 × 10⁻¹ | 4.05 × 10⁻¹ | 1.05 | 12.19 | 13.89 | 12.22 | 8.75 | 6.40 | 3.02 | 1.20 | 0.51 | 0.22 | 0.20 | 0.20 |
| log k1 + 3 | 1.60 | 2.01 | 2.61 | 3.02 | 4.08 | 4.14 | 4.08 | 3.94 | 3.80 | 3.48 | 3.10 | 2.70 | 2.30 | 2.30 | 2.30 |

### Equations (rate law for step 1)

- (8) rate = k1[CA]{[I] + [II] + [III] + [IV]} — I = ⁺H3N/SH, II = ⁺H3N/S⁻, III = H2N/SH, IV = H2N/S⁻ (all COO⁻).
- (9) k1 = {k1′[H⁺]² + (K1[H⁺]/2)(k2′ + k3′) + k4′K1K2} / ([H⁺]² + K1[H⁺] + K1K2), assuming [II] = [III].
- **(10) k1 = (5.50 × 10⁻³ [H⁺]² + 8.40 K1[H⁺] + 0.05 K1K2) / ([H⁺]² + K1[H⁺] + K1K2)**, K1 = 4.07 × 10⁻⁹, K2 = 6.98 × 10⁻¹¹.

So the fitted species constants are k1′ (I) = 5.50 × 10⁻³, (k2′ + k3′)/2 = 8.40 (II/III), k4′ (IV) = 0.05 M⁻¹ s⁻¹ (**F**).
The authors say species IV "was supposed to have greater reactivity" than II/III, so the fall above pH 9.5 is
attributed to a **change of rate-determining step to step 2 (thiazolidine formation)**, not to IV being unreactive —
i.e. k4′ = 0.05 is an empirical fit coefficient, not a microscopic rate constant.

### Table 3. Time course at pH 8.0, 25 °C, I = 0.1 (polarography; "mole" = mol/L)

| t (s) | cinnamaldehyde | cysteine | monoadduct | diadduct |
|---|---|---|---|---|
| 0 | 3.30 × 10⁻⁴ | 6.60 × 10⁻⁴ | 0 | 0 |
| 12 | 3.12 × 10⁻⁴ | 6.24 × 10⁻⁴ | 1.5 × 10⁻⁵ | 0.18 × 10⁻⁴ |
| 24 | 2.95 × 10⁻⁴ | 5.93 × 10⁻⁴ | 3.00 × 10⁻⁶ | 0.32 × 10⁻⁴ |
| 36 | 2.88 × 10⁻⁴ | 5.86 × 10⁻⁴ | 1.00 × 10⁻⁵ | 0.32 × 10⁻⁴ |
| 54 | 2.82 × 10⁻⁴ | 5.60 × 10⁻⁴ | 6.00 × 10⁻⁶ | 0.50 × 10⁻⁴ |
| 140 | 2.52 × 10⁻⁴ | **8.00 × 10⁻⁴** | 0.70 × 10⁻⁴ | (blank) |
| 230 | 2.35 × 10⁻⁴ | 4.76 × 10⁻⁴ | 6.00 × 10⁻⁶ | 0.89 × 10⁻⁴ |
| 350 | 2.13 × 10⁻⁴ | 4.44 × 10⁻⁴ | **1.80 × 10⁻⁴** | 0.99 × 10⁻⁴ |
| 590 | 1.90 × 10⁻⁴ | 4.08 × 10⁻⁴ | **2.80 × 10⁻⁶** | 1.12 × 10⁻⁴ |
| 720 | 1.80 × 10⁻⁴ | 3.92 × 10⁻⁴ | 3.20 × 10⁻⁵ | 1.18 × 10⁻⁴ |
| 900 | 1.70 × 10⁻⁴ | 3.75 × 10⁻⁴ | 3.50 × 10⁻⁵ | 1.25 × 10⁻⁴ |
| 1,200 | 1.67 × 10⁻⁴ | 3.60 × 10⁻⁴ | 2.60 × 10⁻⁵ | 1.37 × 10⁻⁴ |
| 1,500 | 1.60 × 10⁻⁴ | 3.45 × 10⁻⁴ | 2.50 × 10⁻⁵ | 1.45 × 10⁻⁴ |

### Table 4. Rate and equilibrium constants, pH 8.0, I = 0.1

| constant (units as printed) | 13.0 °C | 25.0 °C | 32.0 °C | 38.0 °C |
|---|---|---|---|---|
| k1 (mole⁻¹ sec⁻¹) | 5.852 | 12.19 | 17.45 | 25.08 |
| k−1 (sec⁻¹) | 35.34 | 116.1 | 194.4 | 319.8 |
| k2 (mole⁻¹ sec⁻¹) | 78.33 | 122.3 | 148.4 | 184.9 |
| "k−1" (mole⁻¹ sec⁻¹) [= k−2; see flags] | 5.820 × 10⁻⁵ | 8.185 × 10⁻⁵ | 9.855 × 10⁻⁵ | 11.49 × 10⁻⁵ |
| K1 | 16.56 × 10⁻² | 10.50 × 10⁻² | **7.840 × 10⁻²** | 7.840 × 10⁻² |
| K2 | 13.56 × 10⁵ | 14.94 × 10⁵ | 15.06 × 10⁵ | 16.09 × 10⁵ |
| K1·K2 | 2.229 × 10⁵ | 1.569 × 10⁵ | 1.352 × 10⁵ | 1.262 × 10⁵ |

### Table 5. Activation and thermodynamic parameters, pH 8.0, 25 °C, I = 0.1 (energies kcal/mol, entropies e.u.)

| | step (1) | step (−1) | step (2) | step (−2) | overall |
|---|---|---|---|---|---|
| Ea | 9.17 | 14.9 | 6.11 | 5.96 | |
| ΔH‡ | 8.58 | 14.3 | 5.52 | 5.37 | |
| ΔS‡ | −24.9 | −1.16 | −30.5 | −40.4 | |
| ΔG‡ | 16.0 | 14.65 | 14.6 | 17.4 | |
| ΔH | −5.39 (step 1) | | 1.18 (step 2) | | −4.21 |
| ΔS | −22.4 | | 32.2 | | 9.8 |
| ΔG | 1.31 | | −8.42 | | −7.11 |

Proposed step-2 mechanism (scheme on p. 28): the monoadduct aldehyde + cysteinate (H2N, S⁻) → (Schiff base / hemithioacetal
intermediate, "at least one intermediate in the second step" because ΔG + ΔG‡−2 > ΔG‡1) → thiazolidine + H2O. Both steps
"partially rate determining at pH 8.0". No pH dependence of k2 is measured.

## 4. Arithmetic (mine)

1. **Internal consistency of Table 4.** K1 = k1/k−1: 5.852/35.34 = 0.1656 ✔ (13 °C); 12.19/116.1 = 0.1050 ✔; 25.08/319.8
   = 0.0784 ✔ (38 °C); **17.45/194.4 = 0.0898 at 32 °C, printed 0.0784 (copied from 38 °C)**; K1K2 at 32 °C printed
   1.352 × 10⁵ = 0.0898 × 1.506 × 10⁶ ✔, so the true K1(32 °C) = 8.98 × 10⁻². K2 = k2/k−2: 122.3/8.185 × 10⁻⁵ = 1.494 × 10⁶ ✔
   — this holds only if k−2 is in **s⁻¹**, not the printed "mole⁻¹ sec⁻¹". K1·K2 units M⁻² (two cysteines).
   Overall ΔG = −7.11 kcal/mol → K = exp(7110/(1.987 × 298.15)) = 1.6 × 10⁵ ✔ K1K2.
2. **Arrhenius from Table 4 (13 → 38 °C, Δ(1/T) = 2.808 × 10⁻⁴ K⁻¹):** k2: ln(184.9/78.33) = 0.859 → **Ea ≈ 6.1 kcal/mol
   (25 kJ/mol)** ✔ Table 5. k−2: ln(11.49/5.820) = 0.680 → 4.8 kcal/mol (printed 5.96). van 't Hoff of K2: ln(16.09/13.56)
   = 0.171 → **ΔH2 ≈ +1.2 kcal/mol** ✔ printed 1.18, but **inconsistent with Ea2 − Ea−2 = 0.15** (printed). k1: ln(25.08/5.852)
   = 1.455 → 10.3 kcal/mol (printed 9.17).
3. **Mass balance of Table 3:** [CA]0 = [CA] + mono + di. At 350 s, 3.30 − 2.13 − 0.99 = 0.18 × 10⁻⁴ → mono = **1.8 × 10⁻⁵
   (printed 1.80 × 10⁻⁴ is a typo)**; at 590 s, 3.30 − 1.90 − 1.12 = 0.28 × 10⁻⁴ → **2.8 × 10⁻⁵ (printed 2.80 × 10⁻⁶)**; at
   1500 s CA and Cys balances both close (3.30 ✔; 3.45 + 0.25 + 2 × 1.45 = 6.60 ✔). **Cysteine 8.00 × 10⁻⁴ at 140 s exceeds the
   charge** (likely 5.0 × 10⁻⁴); monoadduct values at 24-230 s scatter by 10× — treat Table 3 as noisy.
4. **Two methods disagree for step 1 at pH 8-9:** k1 = 3.52 (spectro) vs 12.19 (polarography) at pH 8.0; 3.79 vs 13.89 at
   pH 8.5 — ×3.5. Below pH 7.5 they agree within 30 %. The Table 4 rates are the polarographic set.
5. **Timescale for step 2 (transfer test):** at [Cys] = 1 mmol/L, k2[Cys] = 0.12 s⁻¹ (t½ ≈ 6 s) at pH 8, 25 °C; k−2 =
   8 × 10⁻⁵ s⁻¹ (t½ ≈ 2.4 h). Compare Kallen formaldehyde at pH 6.2: K_app ~1.5 × 10⁶ M⁻¹ and hydrolysis t½ ~0.4 h [Z].
   **Same order for K (10⁶ M⁻¹) for two non-sugar aldehydes**; an aldose's open-chain fraction would scale k_f and K down,
   not k_r.

## 5. What the repo could take

| item | value | conditions | grade |
|---|---|---|---|
| k2, aliphatic aldehyde + Cys → thiazolidine | 122 M⁻¹ s⁻¹ (78-185 over 13-38 °C) | pH 8.0, I 0.1 | **F** (fit of eq 7 to Table 3-type data) |
| k−2, thiazolidine ring-opening/hydrolysis | 8.2 × 10⁻⁵ s⁻¹ (5.8-11.5 × 10⁻⁵) | pH 8.0, I 0.1 | **F** (units corrected, !) |
| K2 | 1.49 × 10⁶ M⁻¹; ΔH ≈ +1.2 kcal/mol | pH 8.0, 13-38 °C | **F** |
| Ea(k2) ≈ 6.1 kcal/mol; Ea(k−2) 4.8-6.0 kcal/mol | | pH 8.0, 13-38 °C, 4 points | **F/Z** |
| step-1 (Michael) pH-rate profile, eq 10 | k1′ 5.5 × 10⁻³, (k2′+k3′)/2 8.4, k4′ 0.05 M⁻¹ s⁻¹ | 25 °C, I 0.1 | **F** — not for TTCA |
| k1 obs vs pH (Tables 1-2) | | 25 °C | **M** — not for TTCA |
| **directional/bound:** thiazolidine formation from a free aliphatic aldehyde at pH 8 is seconds-fast at mM cysteine and nearly temperature-insensitive (Ea ~6 kcal/mol); a fast-equilibrium treatment of TTCA at 90-140 °C is consistent | | | **Z** |
| ! k2 measured at ONE pH (8.0); no microscopic species assigned for step 2 beyond the proposed cysteinate attack | | | **!** |
| ! Table 4 row 4 labelled "k−1" with units M⁻¹ s⁻¹; it is k−2 in s⁻¹ (K2 = k2/k−2 only then) | | | **!** |
| ! K1(32 °C) misprinted (0.0784 → 0.0898) | | | **!** |
| ! Table 3 typos (350 s, 590 s monoadduct; 140 s cysteine) | | | **!** |
| ! ΔH2 (+1.18) vs Ea2 − Ea−2 (0.15) inconsistent | | | **!** |
| ! Cysteine pKas from Wrathall 1964, which Kallen 1971 questions | | | **!** |
| ! Monoadduct is a 3-thio-substituted aldehyde, not an aldose; no ring-chain tautomerism | | | **!** |

## Use in the model

- **B59b rates:** if the engine wants a k_f for Cys + PENT → TTCA with a declared band, Kim's k2 (≈ 10² M⁻¹ s⁻¹, pH 8, 25 °C)
  is the corpus's only measured-rate analogue for an aliphatic aldehyde; scaled by an aldose open-chain fraction it would
  still give sub-minute equilibration at 0.1-0.4 M and 90 °C. Use as a **bound supporting the fast-equilibrium assumption**,
  not as a constant. Ea ≈ 25 kJ/mol is a weak temperature dependence; do not carry it to TTCA without a band.
- **pH:** this paper does not give the pH dependence of thiazolidine formation; use Kallen 1971 for the shape.

# Kabyemela, Adschiri, Malaluan & Arai 1999 — EXTRACTION
### *"Glucose and Fructose Decomposition in Subcritical and Supercritical Water: Detailed Reaction Pathway, Mechanisms, and Kinetics"*
### Flow reactor, pure water, 300–400 °C, 25–40 MPa, 0.02–2 s. Glucose and fructose break by retro-aldol into erythrose + glycolaldehyde (C4 + C2) and into 2 glyceraldehyde (C3 + C3). Thirteen first-order constants at four conditions; **no activation energy and no pre-exponential is printed**. Read 2026-10-27.

**Full read of `data/articles/Kabyemela1999.pdf`** (text layer, 8 pages; Figures 6 and 9 and Table 3
checked on 90-dpi rasters). The earlier papers it builds on (Kabyemela 1997 *IECR* 36:1552, 2025,
5063) are **not on disk**.

**Provenance:** **[M]** measured · **[F]** fitted by the authors · **[Q]** qualitative ·
**[D]** derived by this extraction from printed numbers · **[Z]** extrapolated by this extraction
outside the paper's range · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Bernard M. Kabyemela*, Tadafumi Adschiri, Roberto M. Malaluan, Kunio Arai |
| Venue | *Ind. Eng. Chem. Res.* 1999, 38 (8), 2888–2895 |
| DOI | 10.1021/ie9806390 (printed p. 2888 footer) |
| Affiliation | Department of Chemical Engineering, Tohoku University, Sendai |
| Dates | received 6 Oct 1998; revised 22 Apr 1999; accepted 30 Apr 1999; web 17 Jul 1999 |
| PDF character | 8 pp.; Tables 1–3; Figures 1–9 (Fig. 6 mechanism, Fig. 9 kinetic scheme; Figs. 4, 7, 8 data, not digitised) |

---

## 1. THE ONE-PARAGRAPH ANSWER

In sub- and supercritical water (300–400 °C, 25–40 MPa, residence 0.02–2 s, no catalyst, no buffer,
no amine), glucose isomerises to fructose, dehydrates to 1,6-anhydroglucose, and **cleaves C–C by
retro-aldol**: the open-chain aldose (or the 2,3-enediol) splits C2–C3 into **erythrose +
glycolaldehyde**, and the ketose (or the 1,2-enediol) splits C3–C4 into **2 glyceraldehyde**, which
isomerises to dihydroxyacetone; both dehydrate to **pyruvaldehyde (methylglyoxal)**. Fructose does
the same except it makes no 1,6-anhydroglucose and barely isomerises back. Table 3 prints 13 rate
constants (s⁻¹) at 300 °C/25 MPa, 350 °C/25 MPa, 400 °C/30 MPa and 400 °C/40 MPa. At 300 °C the
glucose C–C cleavage constants are k_GE = 0.21 s⁻¹ (→ erythrose + glycolaldehyde) and k_GGLY = 0.05
s⁻¹ (→ 2 glyceraldehyde), together **55 %** of glucose loss [D]. Glycolaldehyde was identified by
¹H-NMR but **not quantified**; k_GE is fitted on erythrose. **No Arrhenius parameters are printed**:
everything at 145 °C below is this extraction's two- or three-point extrapolation across 155 K.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Reactor | continuous flow, stainless steel, ID 0.077 and 0.118 cm, in a KNO3 + KNO2 salt bath; residence time varied by reactor length |
| Feed | aqueous glucose ≈ 0.6 wt %, or fructose or 1,6-anhydroglucose ≈ 0.25 wt %, at ≈ 5 mL/min, mixed with preheated water at 20 mL/min just before the reactor (so ≈ 5× dilution; ≈ 6.7 mM glucose in-reactor [D], consistent with the ≈ 6.5–7 mM starting points in Fig. 8 [Q]) |
| Quench | cool water injected at 14 mL/min at the exit plus a cooling jacket, < 60 °C |
| Conditions | 300 and 350 °C at 25 MPa (subcritical); 400 °C at 30 and 40 MPa (supercritical); 0.02–2 s |
| pH / buffer | **none stated; pure water** [!] |
| Chemicals | Wako: glucose 99+ %, fructose 99+ %, glyceraldehyde 97+ %, erythrose 60 %, pyruvaldehyde 40 %, dihydroxyacetone 99+ %; TCI: 1,6-anhydroglucose 99+ % |
| Analysis | HPLC: Shodex KS 802, water 1 mL/min, 353 K, RI (quantitative) + UV 290 nm (C=O confirmation); ¹H NMR (Asahi Chemical) to confirm fructose, erythrose, glyceraldehyde, DHA, pyruvaldehyde, glycolaldehyde, 1,6-anhydroglucose; TOC |
| Carbon balance | 95–100 % (glucose runs), 80–95 % (fructose runs); > 90 % at short residence, falling at longer times as acetic/formic acids (detected, not quantified) accumulate; gas negligible |
| Fitting | MLAB, Gear third-order integrator; sequential: fructose constants first (k_FTOT from −ln(1 − X) vs t, Fig. 7, first order), then glucose constants with fructose and glyceraldehyde/DHA constants fixed; r = 0.90–0.98; standard deviations "using the procedure explained by S. W. Benson" |

---

## 3. THE KINETIC SCHEME (Figure 9, eqs 1–10, as printed)

| symbol | reaction | stoichiometry in the ODEs |
|---|---|---|
| k_GF | glucose → fructose | eq 8 |
| k_GA | glucose → 1,6-anhydroglucose | eq 5 |
| **k_GE** | **glucose → erythrose (+ glycolaldehyde)** | eq 6: d[e]/dt = k_GE[g] + k_FE[f] − k_E[e]; glycolaldehyde not in the ODEs |
| **k_GGLY** | **glucose → 2 glyceraldehyde** | eq 7: d[gly]/dt = 2k_GGLY[g] + 2k_FGLY[f] − … |
| k_FE | fructose → erythrose (+ glycolaldehyde) | eqs 1, 6 |
| k_FGLY | fructose → 2 glyceraldehyde | eqs 2, 7 |
| k_FACID | fructose → acetic + formic acids | eq 4; obtained as k_FTOT − (k_FE + k_FGLY) |
| k_E | erythrose → liquid products (acids) | eqs 1, 6 |
| k_A | 1,6-anhydroglucose → liquid products (acids) | from separate 1,6-anhydroglucose runs |
| k_GLYD / k_DGLY | glyceraldehyde ⇄ dihydroxyacetone | from Kabyemela 1997 (ref. 9) |
| k_GLYP / k_DP | glyceraldehyde / DHA → pyruvaldehyde (MGO) | from Kabyemela 1997 (ref. 9) |

Mechanism (Figure 6, [Q]): open-chain glucose **2** → retro-aldol at C2–C3 → erythrose + glycolaldehyde;
1,2-enediol **3** → C3–C4 cleavage ("double bond rule", β to the C=C) → glyceraldehyde; fructose keto
form **4** → retro-aldol → glyceraldehyde; 2,3-enediol **7** → C2–C3 cleavage → erythrose +
glycolaldehyde. DHA is a **secondary** product (no intercept in selectivity vs conversion, Fig. 5,
350 °C) [Q]. 5-HMF "possible … not significant" at these residence times [Q]. Erythrose and
1,6-anhydroglucose give only acetic and formic acids, no glyceraldehyde or DHA (Table 2) [M/Q].

---

## 4. EVERY PRINTED NUMBER

### 4.1 Table 3 — rate constants, s⁻¹, ± "standard deviation" [F] (p. 2892)

| k (s⁻¹) | 300 °C, 25 MPa | 350 °C, 25 MPa | 400 °C, 30 MPa | 400 °C, 40 MPa |
|---|---|---|---|---|
| k_GA | 0.010 ± 0.001 | 0.020 ± 0.001 | 0.080 ± 0.007 | 0.500 ± 0.045 |
| k_E | 0.100 ± 0.004 | 0.550 ± 0.022 | 5.000 ± 0.450 | 4.000 ± 0.360 |
| k_A | 0.010 ± 0.001 | 0.040 ± 0.002 | 0.310 ± 0.028 | 1.440 ± 0.130 |
| k_FGLY | 0.100 ± 0.004 | 0.600 ± 0.024 | 6.500 ± 0.585 | 5.000 ± 0.45 |
| k_FACID | 0.180 ± 0.007 | 0.700 ± 0.028 | 10.400 ± 0.936 | 4.800 ± 0.432 |
| k_FE | 0.100 ± 0.004 | 0.800 ± 0.032 | 8.000 ± 0.720 | 6.000 ± 0.540 |
| **k_GE** | **0.210 ± 0.008** | **0.950 ± 0.038** | 18.100 ± 1.629 | 12.000 ± 1.080 |
| **k_GGLY** | **0.050 ± 0.002** | **0.200 ± 0.008** | 1.000 ± 0.090 | 2.500 ± 0.225 |
| k_GF | 0.200 ± 0.008 | 0.640 ± 0.026 | 7.000 ± 0.630 | 6.000 ± 0.540 |
| k_GLYD ᵃ | 0.400 ± 0.016 | 1.380 ± 0.052 | 7.150 ± 0.644 | 15.700 ± 1.413 |
| **k_GLYP** ᵃ | 0.190 ± 0.008 | 0.940 ± 0.038 | 4.600 ± 0.414 | 7.420 ± 0.668 |
| k_DGLY ᵃ | 0.030 ± 0.001 | 0.200 ± 0.008 | 1.040 ± 0.094 | 2.450 ± 0.794 |
| k_DP ᵃ | 0.170 ± 0.007 | 0.560 ± 0.022 | 1.200 ± 0.108 | 2.550 ± 0.230 |

ᵃ "Reaction rate from previous work (Kabyemela et al. 9)".

### 4.2 Other printed numbers

| quantity | value | location |
|---|---|---|
| Glyceraldehyde → DHA vs reverse | "always about 7 times larger" | p. 2892 [Q] |
| Glyceraldehyde → pyruvaldehyde vs DHA → pyruvaldehyde | glyceraldehyde "always greater" | p. 2892 [Q] |
| Pressure effect | subcritical → supercritical "did not result in significant changes in the rate constant" (Arrhenius plot of glucose decomposition, from ref. 8, not reproduced) | p. 2894 [Q] |
| Erythrose stability | "relatively stable, even at higher temperatures, notably in the supercritical region" | p. 2893 [Q] |
| Fructose decomposition order | −ln(1 − X) vs t linear at all four conditions: first order | p. 2893, Fig. 7 [M, figure] |

### 4.3 Derived from Table 3 [D]

Share of glucose loss going to C–C cleavage, (k_GE + k_GGLY)/(k_GA + k_GE + k_GF + k_GGLY):
**0.55** at 300 °C, 0.64 at 350 °C, 0.73 at 400 °C/30 MPa, 0.69 at 400 °C/40 MPa. The erythrose +
glycolaldehyde channel alone: 0.45, 0.52, 0.69, 0.57.

Two-point Arrhenius on the two subcritical conditions (300 and 350 °C, both 25 MPa), and three-point
adding 400 °C/30 MPa:

| k | Ea, 2-pt (kJ/mol) | A, 2-pt (s⁻¹) | Ea, 3-pt (kJ/mol) |
|---|---|---|---|
| k_GE (Glc → erythrose + glycolaldehyde) | 89.6 | 3.1e7 | 141.4 |
| k_GGLY (Glc → 2 glyceraldehyde) | 82.3 | 1.6e6 | 95.7 |
| k_FE (Fru → erythrose + glycolaldehyde) | 123.5 | 1.8e10 | 140.1 |
| k_FGLY (Fru → 2 glyceraldehyde) | 106.4 | 5.0e8 | 133.1 |
| k_GF (Glc → Fru) | 69.1 | 4.0e5 | 112.8 |
| k_GA (Glc → 1,6-anhydroglucose) | 41.2 | 56.5 | 66.0 |
| k_GLYP (glyceraldehyde → MGO) | 95.0 | 8.6e7 | 102.0 |
| total glucose loss | 80.1 | 9.3e6 | — |

The 2-point Ea are low for sugar chemistry (41–124 kJ/mol against 113–159 in Brands 2002c at 90–130
°C) and the 3-point values differ from them by up to 52 kJ/mol: the 400 °C points are supercritical,
where density and solvent properties change, so a single Arrhenius line through 300–400 °C is not
physically one regime.

---

## 5. FLAGS [!]

1. **No Arrhenius parameters printed.** The abstract and conclusion promise kinetics "to predict the
   product distribution … at different conditions", but only four-condition rate tables exist; the
   Arrhenius analysis cited (ref. 8) covers total glucose decomposition only and is not reproduced.
2. **The ± is not a fit uncertainty.** Every entry at 300/350 °C carries ≈ 4 % and every 400 °C entry
   ≈ 9 % (e.g. 0.100 ± 0.004, 5.000 ± 0.450, 18.100 ± 1.629): a fixed fractional error from the Benson
   procedure, not a confidence interval from the fit. One exception, k_DGLY at 400 °C/40 MPa (2.450 ±
   0.794, 32 %), is probably a typo.
3. **Pressure is not innocuous** in the table even though the text says it is: at 400 °C, k_GA is
   6.3× higher at 40 MPa than at 30 MPa, k_GGLY 2.5×, k_FACID 0.46×.
4. Glycolaldehyde is **identified, not quantified**; the C2 rate is inferred from erythrose (k_GE, k_FE)
   on the assumption of 1:1 co-formation. Erythrose is itself consumed (k_E), so k_GE depends on k_E,
   which was "simulated … until it gave a reasonable value".
5. pH is never stated. Unbuffered water at 300–400 °C has a very different ion product and dielectric
   constant from a pH 4–7 food buffer at 100–150 °C; acid/base catalysis is folded into the constants.
6. Text typos: eq 7 writes "k_FGLYA"; Fig. 9 labels fructose → acids "k_FA" (k_FACID in Table 3 and
   notation); notation lists "GTO" for overall glucose decomposition, unused; p. 2894 lists k_GLYP
   twice ("k_GLYD, k_GLYP, k_DGLY, and k_GLYP"); the glyceraldehyde/DHA constants are credited to ref.
   10 (an AIChE 1996 presentation) on p. 2892 and to ref. 9 in the Table 3 footnote.
7. "About 7 times" (k_GLYD/k_DGLY) holds at 350 and 400/30 (6.9) but not at 300 °C (13.3) or
   400/40 (6.4) [D].
8. Abstract says 25–40 MPa for all temperatures; the objectives paragraph says "300–350 °C and 25–40
   MPa, 400 °C at 30–40 MPa"; Table 3 has only 25 MPa below 400 °C.

---

## 6. WHAT THE MODEL CAN TAKE FROM IT

**The need:** `k_glc_ha`, the sulfur lane's `r_glc_c2c3`, Glc → HA + MGO + FRAG_C, first order,
1/min, reference 145 °C (418.15 K). Frozen B9: log10 k = −8.59 (2.6e-9 /min); B47 refit: −4.95
(1.1e-5 /min); unidentified.

- **This is the only paper on disk so far with a rate constant for glucose → glycolaldehyde.**
  But the product split is **C4 + C2 (erythrose + glycolaldehyde, k_GE)** or **C3 + C3 (2
  glyceraldehyde → DHA → MGO, k_GGLY)**, never C2 + C3 + C1 in one event. `r_glc_c2c3` is a
  stoichiometric composite of the two measured channels. If `k_glc_ha` is read as the **HA-producing**
  rate, k_GE is the analogue (1 HA per glucose; its partner is erythrose, which goes to acids, not
  MGO). If it is read as the **MGO-producing** rate, k_GGLY is (2 C3 per glucose, then k_GLYP/k_DP
  to MGO). No amine present; no buffer; pure water.
- **Units:** Table 3 is s⁻¹; × 60 for the model's 1/min.
- **Extrapolation to 145 °C [Z]** (155 K below the lowest measured temperature; the 1/T lever is 4.6×
  the measured 300–350 °C span):

  | channel | 2-pt (300/350 °C) at 145 °C | 3-pt (+400 °C/30 MPa) at 145 °C |
  |---|---|---|
  | k_GE, Glc → erythrose + **HA** | 1.2e-2 /min (log10 −1.93) | 1.6e-4 /min (log10 −3.80) |
  | k_GGLY, Glc → 2 glyceraldehyde (→ MGO) | 5.0e-3 /min (log10 −2.30) | 1.6e-3 /min (log10 −2.79) |
  | k_FE, Fru → erythrose + HA | 4.0e-4 /min | 1.0e-4 /min |
  | k_FGLY, Fru → 2 glyceraldehyde | 1.5e-3 /min | 1.7e-4 /min |

  The two fits disagree by up to **1.9 decades** for k_GE. **Caveat:** the 2-point Ea values are too
  low to be barrier-like (k_GA's 41 kJ/mol is not plausible for a dehydration); the solvent's ion
  product, density and dielectric constant all change between 145 °C and 300–400 °C, and these
  changes are absorbed into the apparent Ea. A plausibility check points the same way: the 2-point
  extrapolation puts total glucose loss in pure water at 120 °C near 0.03 /min, about 2.5× the
  0.012 /min Brands 2002c measures in **phosphate**-catalysed pH 6.8 buffer at 120 °C (Table 6:
  k1 + k3 + k14), whereas unbuffered water should be slower. So the 2-point numbers are high-side.
- **Verdict: a bound, not an anchor.** No printed number is at food temperature, no Ea is printed, the
  HA channel is never quantified directly, and the derived temperature dependence spans 1.9 decades
  at 145 °C depending on which points are used. What survives: (1) the **mechanistic claim** [Q]
  that hexose C–C cleavage is retro-aldol from the open chain and the enediols, giving C2 + C4 and
  C3 + C3, which is the structural basis for a hexose → HA step at all; (2) a **soft ceiling**:
  amine-free hexose fragmentation to HA at 145 °C is at most ~10^−2 /min, and plausibly
  10^−4–10^−3 /min. Both B9 (−8.59) and B47 (−4.95) are below even the 3-point value (−3.80), so
  this ceiling does not bind either fit. It also cannot discriminate between them. As with Brands
  2002c, amine catalysis in a cysteine pot could raise the real rate above a sugar-only value, so the
  ceiling is soft.

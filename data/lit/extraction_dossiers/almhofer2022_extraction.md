# Almhofer, Bischof, Madera & Paulik 2022 — EXTRACTION (furfural degradation in water, 125-200 °C, uncatalysed and in 0.1 M of ten acids; overall Ea 49.8 kJ/mol, order 0.8 to 1.3)

**Source on disk:** `data/articles/almhofer2022.pdf` (the publisher's open-access PDF, 17 pp.). Pages 3 and
5-10 and 13 (every page carrying a number used below) read and checked by eye on 2026-10-09 from the page
images, cross-checked against `pdftotext -layout`; the rest read from the text layer. Figure 2 (p. 6) and Figure 3 (p. 7) were re-read from 250-300 dpi crops. The
Supporting Information (Sections S1-S4, Tables S3, Figure S2) is **not on disk**. Written for the
`k_fur_decay` comparison.

| field | value |
|---|---|
| Title | "Kinetic and mechanistic aspects of furfural degradation in biorefineries" |
| Authors | Lukas Almhofer, Robert H. Bischof, Martin Madera, Christian Paulik (Wood K plus / JKU Linz / Lenzing AG) |
| Venue | Can. J. Chem. Eng. 101:2033 (volume and first page as cited in the request; the PDF on disk is the early-view version, header "Can J Chem Eng. 2022;1-17", received 5 Jan, accepted 14 April 2022; final issue pagination not printed in it) |
| DOI | 10.1002/cjce.24593 |

## 1. Methods

- **Reactors (p. 3).** Stainless-steel batch reactors, ~8 mL, 7 mL charge, stirred 500 rpm in a heated
  aluminium block; quenched in ice. At least duplicate. **The metal surface is in contact with the
  solution**; the authors name metal-surface catalysis as one possible cause of the sub-first order at
  150 °C (p. 5) and H2SO4 attacking the steel walls (p. 8).
- **Uncatalysed arm (Section 3.1, p. 5-7):** furfural in water, 5, 10, 25, 50 g/L (**52, 104, 260,
  520 mmol/L, derived here**, MW 96.08), 125, 150, 175, 200 °C, to 360 min. pH not printed.
- **Acid arm (Section 3.2, p. 7-8):** 0.1 mol/L of each acid, 1 mass % furfural, 150 °C. pH at 150 °C
  calculated from Ka(T) (Eq. 3, Table 1), not measured hot.
- **Analysis.** Furfural by HPLC-UV 277 nm (external calibration); formate by ion chromatography.
- **Kinetics.** Order from initial rates (conversion < 32 %, R² > 0.98). Overall Ea assumes **first
  order**, at 10 g/L. Mechanistic model (Eq. 4-6) fitted in MATLAB at 150 °C only.

## 2. Findings that matter

**Uncatalysed furfural loss, Arrhenius (Figure 2B inset, p. 6; text p. 7).** Inset printed:
intercept **5.74E+00 ± 9.95E-01**, slope **-5.99E+03 ± 4.31E+02** (K), R² **0.990**. Text: "overall
activation energy of **49.8 kJ mol⁻¹**", "similar range of published values without an acid catalyst,
ranging from 44.2 to 58.8"; acid-catalysed literature 48.1-110.3 (cited, ref. 39).

Units of k are **not printed**. Time axis is min (Fig. 2A); the nomenclature (p. 15) gives t in s. The
inset intercept read as ln(k / min⁻¹) reproduces Figure 1A's directly plotted conversion (150 °C, 10 g/L,
X ≈ 7 % at 360 min gives -ln(0.93)/360 = 2.0e-4 min⁻¹ against 2.2e-4 from the fit), so **k is in min⁻¹
(inferred here, not printed)**. Figure 2A's y-axis is labelled ln[1-X(FU)] but its values (to -40) are
consistent only with 100 × ln(1-X); a labelling slip in the figure.

| T (°C) | ln k, read from Fig. 2B (approx.) | k from the printed fit, min⁻¹ (derived here) |
|---|---|---|
| 125 | -9.21 | 9.1e-5 |
| 150 | -8.52 | 2.2e-4 |
| 175 | -7.67 | 4.9e-4 |
| 200 | -6.82 | 9.9e-4 |
| 145 (interpolated) | — | **1.87e-4** (log10 -3.73) |

Derived here: Ea = 5990 × 8.314 = 49.8 ± 3.6 kJ/mol; A = e^5.74 = **311 min⁻¹** (e^(5.74 ± 0.995):
115-841 min⁻¹).

**Reaction order (p. 5, Fig. 1C/D):** initial-rate orders **0.8 at 150 °C and 1.3 at 200 °C**. At 150 °C
conversion falls with rising initial concentration (Fig. 1A: 5 g/L ~11 %, 10 g/L ~7 %, 25 and 50 g/L
~6 % at 360 min, read from graph, approx.); at 200 °C it rises (Fig. 1B, 240 min: 5 g/L ~16.5 %, 10 g/L
~21 %, 25 g/L ~27.5 %, 50 g/L ~31 %, read from graph, approx.).

**Acid arm, 150 °C, 1 mass % furfural, conversion at 360 min (Figure 3, p. 7; read from graph, approx.):**

| acid (0.1 M) | pH at 150 °C (printed) | X(FU) at 360 min |
|---|---|---|
| H2SO4 | 1.00 | ~33 % |
| LSA | 1.48 (room-T value, flagged) | ~44.5 % |
| oxalic | 1.56 | ~17 % |
| H3PO4 | 1.97 | ~12 % |
| citric | 2.19 | ~11 % |
| SO2 | 2.22 | ~23 % |
| formic | 2.61 | ~7 % |
| succinic | 2.77 | ~8.5 % |
| acetic | 3.09 | ~7.6 % (last point ~370 min) |
| propionic | 3.15 | ~6 % |
| Na2H-citrate | 6.32 | ~11.5 % |

Conversion rises monotonically with falling pH except SO2 and LSA (bisulfite adducts, condensation) and
formic acid (reversible formic-acid pathway). At pH 2.5, 150 °C, 10 g/L, acetic/succinic/phosphoric acids
gave 6.23, 6.29, 6.35 (p. 8; % implied; the time point is not printed). Formic-acid selectivity 8-35 %
(p. 9). Formic-acid pathway Ea **27.0 kJ/mol**, first order assumed (p. 10; Arrhenius plot in SI).

**Mechanistic model at 150 °C (Eq. 6, p. 13), as printed:**
d[FU]/dt = -2.2 × 10⁻³ [FU]² - 2.6 × 10⁻³ [H⁺]^0.58 [FU] + 7.3 × 10⁻² [H⁺]^0.58 [FA].
k1 (second-order "direct polymerisation", uncatalysed) = 2.2e-3; k2 = 2.6e-3; k2′ = 7.3e-2; h = 0.58.
**Units of none of these are printed** (concentration in mol/L per nomenclature; time base ambiguous,
min or s). SO2, H2SO4 and LSA were excluded from the fit.

## 3. What it means for the model

The engine's coordinate is `k_fur_decay` (`r_fur_decay`: FUR -> 5 FRAG_C, first order, 1/min, the "large
unidentified sink", `src/kinetic_core/parameters_sulfur.py` line 2108). Live values from
`results/validation/core_prediction_uncertainty.json`:

| key | centre | distribution / band | reason |
|---|---|---|---|
| `b8.k_fur_decay.log10_k_ref_145C` | **0.470** (k = 2.95 min⁻¹ at 145 °C, derived here) | normal_log10, σ 0.437, band [-10, 0.5] | laplace_covariance_at_b8_optimum; **bound_limited** in data_wishlist §1 (centre 0.03 below the 0.5 ceiling) |
| `b8.decay_Ea_kJ_mol.carbonyl_sink` (shared by `k_fur_decay`, `k_osone_decay`, `k_nf_decay`) | **174.9 kJ/mol** | uniform_band [126.9, 223.0] | unidentified_in_the_fit: declared band capped by the prefactor prior |

**Comparison (derived here).** At 145 °C the engine's furfural sink is 2.95 min⁻¹; Almhofer's water-only
loss at 104 mM is 1.87e-4 min⁻¹. The engine is **4.2 decades faster**. The barrier is 175 kJ/mol against
49.8 measured; Almhofer's value lies **below the engine's sampled band** [126.9, 223.0] altogether.
Extrapolated outside Almhofer's range to 100 °C: 3.3e-5 min⁻¹ (Almhofer) against 6.8e-3 min⁻¹ (engine,
with its 174.9 barrier).

**It does not falsify `k_fur_decay`, and cannot pin it.** The engine's lump is all unidentified furfural
consumption in a Maillard pot (amines, H2S, cysteine, other carbonyls), while Almhofer measures furfural in
water alone, at 50-500 mM, in steel. Water-only thermolysis is a **floor** under that lump: the engine
must stay above it, and it does by four decades. Concentration does not close the gap: with order 0.8 at
150 °C, going from 104 mM to ~1 mM raises the apparent first-order constant by about 104^0.2 = 2.5× (derived
here), not 10⁴. What it does say: (i) the furfural sink in the engine is not furfural self-degradation, so
it must be furfural reacting with something in the pot, and the wishlist's proposed fed-furfural
measurement should be run WITH the pot's co-reactants, not alone, or it will measure this floor;
(ii) a 175 kJ/mol barrier on that lump has no support here; water-only loss runs at 50 kJ/mol, and the
acid-catalysed literature the authors cite tops out at 110 kJ/mol.

## What it does not give

- The units of k anywhere (min⁻¹ inferred, not printed); the per-temperature k values as numbers (only
  the Arrhenius inset and the plot); the pH of the uncatalysed arm.
- Any rate below 125 °C, any rate at Maillard-relevant concentrations (≤ a few mM), any rate with an
  amine, sulfide or thiol present.
- Fit uncertainties on the Eq. 6 parameters; the Section S1-S4 SI (initial-rate calculation, formic-acid
  Arrhenius plot, selectivity model).
- Any check that the steel reactor walls do not contribute (the authors raise it themselves).

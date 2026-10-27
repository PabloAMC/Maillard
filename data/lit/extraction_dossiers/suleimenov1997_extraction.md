# Suleimenov & Seward 1997 — EXTRACTION
### *"A spectrophotometric study of hydrogen sulphide ionisation in aqueous solutions to 350°C"*
### The first ionisation constant of H2S (molal, I → 0) from 25 to 350 °C by UV spectrophotometry. pK1 falls from 6.99 at 25 °C to a flat minimum of 6.49 at 100–150 °C, then rises to 8.89 at 350 °C. Eq. 22 gives pK1(T) at s.v.p.; Table 3 gives ΔH and ΔCp of ionisation (100 bar). Read 2026-10-27.

**Full read of `data/articles/suleimenov1997.pdf`** (12 pages, scan with a text layer). The OCR
garbles the exponent in Eq. 22 ("10-aT 2"), so **Eq. 22 and Tables 1–3 were read off a 200-dpi
raster** of pp. 5191–5194, not the text layer.

**Provenance:** **[M]** printed value · **[Z]** derived by me from printed values, never printed ·
**[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | O. M. Suleimenov, T. M. Seward |
| Title | "A spectrophotometric study of hydrogen sulphide ionisation in aqueous solutions to 350°C" |
| Venue | *Geochim. Cosmochim. Acta* **61** (24), **5187–5198** (1997) |
| DOI | 10.1016/S0016-7037(97)00291-3 (the paper prints the PII "S0016-7037(97)00291-3"; the DOI is the standard Elsevier 10.1016 prefix on that PII) |
| Affiliation | Institut für Mineralogie und Petrographie, ETH-Zentrum, Zürich |
| Dates | received 20 March 1997, accepted in revised form 15 July 1997 |

---

## 1. WHAT IS MEASURED, AND HOW [M]

- **Quantity:** K1 for H2S(aq) = H⁺ + HS⁻, molal scale, from the c.t.t.s. UV band of HS⁻
  (250–266 nm, 16 wavelengths). pK1 is fitted by minimising absorbance residuals against a Gibbs-
  energy-minimisation speciation (Eq. 10). The H2S(aq) contribution at those wavelengths is negligible
  (Fig. 6).
- **Solutions (Table 1):** NaOH 2.0–9.6 × 10⁻³, NaH2PO4 7.3–12.6 × 10⁻³ and H2S 0.20–1.02 × 10⁻³
  mol dm⁻³. Phosphate buffer; pH(20 °C) measured 6.09–8.10, calculated 6.08–8.08 (11 solutions).
  Phosphate's second ionisation at temperature from Mesmer & Baes 1974.
- **Ionic strength:** "always less than 0.03" (p. 5190). Activity coefficients by "the standard
  Debye-Hückel equation", γ(H2S) = 1. The reported K1 is therefore effectively the I = 0
  thermodynamic constant.
- **Cell:** gold-lined Ti-alloy optical cell with sapphire windows, ±1 °C. Pressures held just
  above s.v.p. (100 bar to 300 °C, 110 bar at 300 °C in Table 2, 200 bar at 350 °C). Sulphide
  analysed after each run by methylene blue (667 nm), calibrated against a coulometric H2S
  generator.
- **Pressure correction to s.v.p.:** from partial molar volumes (Sretenskaya 1976 at low T;
  SUPCRT92 at high T). Δlog K < 0.03 for 25–300 °C, 0.12 at 350 °C.

---

## 2. TABLE 2 — FIRST IONISATION CONSTANT OF H2S (molal), 25–350 °C [M]

| t (°C) | p (bar) | pK1 (at p) | pK1 (s.v.p.) | uncertainty |
|---|---|---|---|---|
| 25 | 100 | 6.96 | 6.99 | ±0.022 |
| 50 | 100 | 6.66 | 6.68 | ±0.018 |
| 100 | 100 | 6.46 | 6.49 | ±0.013 |
| 150 | 100 | 6.47 | 6.49 | ±0.019 |
| 200 | 100 | 6.70 | 6.73 | ±0.023 |
| 250 | 100 | 7.16 | 7.19 | ±0.076 |
| 300 | 110 | 7.87 | 7.89 | ±0.076 |
| 350 | 200 | 8.77 | 8.89 | ±0.051 |

---

## 3. EQUATION 22, VERBATIM [M] (p. 5192)

Text: it represents the temperature dependence at saturation water vapour pressure and "should be
used for interpolation purposes only".

> lg K1 = 782.43945 + 0.361261 T − 1.6722 × 10⁻⁴ T² − 20565.7315 / T − 142.741722 ln (T)   (22)

"where T is temperature in Kelvin." lg = log10; pK1 = −lg K1. Valid range: 25–350 °C at s.v.p.

### 3.1 Eq. 22 vs Table 2 (s.v.p. column) [Z]

| t (°C) | Eq. 22 pK1 | Table 2 pK1 (s.v.p.) | Eq. 22 − Table 2 | Table 2 uncertainty |
|---|---|---|---|---|
| 25 | 6.978 | 6.99 | −0.012 | ±0.022 |
| 50 | 6.701 | 6.68 | +0.021 | ±0.018 |
| 100 | 6.467 | 6.49 | −0.023 | ±0.013 |
| 150 | 6.499 | 6.49 | +0.009 | ±0.019 |
| 200 | 6.736 | 6.73 | +0.006 | ±0.023 |
| 250 | 7.188 | 7.19 | −0.002 | ±0.076 |
| 300 | 7.891 | 7.89 | +0.001 | ±0.076 |
| 350 | 8.890 | 8.89 | 0.000 | ±0.051 |

Eq. 22 reproduces Table 2 to ≤ 0.023. At 50 and 100 °C the residual slightly exceeds the stated
point uncertainty, so the smooth curve, not the point, is the better estimate there. Interpolated
[Z]: 120 °C 6.454, 130 °C 6.461, **145 °C 6.486**, 160 °C 6.530, 180 °C 6.617. **The minimum is
6.454 at ≈ 117 °C**, and pK1 stays within 6.45–6.57 between 90 and 170 °C.

---

## 4. THERMODYNAMIC FUNCTIONS OF IONISATION [M]

**Table 3** (at **100 bar**, Clarke–Glew fit, Eq. 23). **b2 = ΔCp = −265.3 J mol⁻¹ K⁻¹ at 298.15 K was
imposed** from Barbero et al. 1982's heat capacities, not fitted. Regression coefficients: b0
−133.0347, b1 76.2421, b2 −265.30, b3 196.1261, b4 −348.9833, b5 −122.9204. Selected rows:

| t (°C) | ΔG° (J/mol) | ΔH° (J/mol) | ΔCp° (J mol⁻¹ K⁻¹) |
|---|---|---|---|
| 25 | 39666 | **22733** | **−265.3** |
| 50 | 41346 | 16448 | −240.1 |
| 75 | 43495 | 10592 | −231.3 |
| 100 | 46062 | 4731 | −240.7 |
| 120 | 48407 | −280 | −262.5 |
| 145 | 51718 | −7379 | −309.1 |
| 150 | 52434 | −8954 | −321.1 |
| 175 | 56310 | −17863 | −395.6 |
| 200 | 60739 | −28944 | −495.2 |

ΔH of ionisation **changes sign near 119 °C** (100 bar), the pK1 minimum (Eq. 22 at s.v.p.:
117 °C). **[!]** Table 3 prints ΔH
= −89146 at 275 °C between −77234 (270) and −86860 (280). The sequence is otherwise monotone, so
the 275 °C entry is probably a misprint (about −82 000 expected). It does not affect ≤ 200 °C.

[Z] For comparison, ΔH and ΔCp implied by Eq. 22 (s.v.p.) by numerical differentiation: 25 °C 24.1
kJ/mol, −316; 100 °C 4.35, −246; 145 °C −7.8, −307; 200 °C −29.1, −488. They are consistent with
Table 3. The authors warn that thermodynamic functions at s.v.p. are not presented because the
pressure derivative terms are omitted.

---

## 5. THE MODEL'S IMPLEMENTATION — `parameters_sulfur.h2s_pk1_suleimenov_seward`

Code (working tree, branch `b51-ph-state-everywhere`, **not edited**):
`log_k1 = 782.43945 + 0.361261*t - 1.6722e-4*t*t - 20565.7315/t - 142.741722*math.log(t)`;
`return -log_k1`, with `t` in K.

- **Coefficients, signs and the 10⁻⁴ exponent: identical to Eq. 22.** Natural log is used for "ln",
  and pK = −lg K. **Implements Eq. 22 exactly.** [Z] Evaluated in the `maillard` env at 25/50/100/
  145/150/200 °C, it equals my independent evaluation of Eq. 22 to machine precision (difference
  0.0).
- **[!] Docstring numbers are slightly off (code is right):** the docstring says "6.97 at 25 C,
  6.45 at 100 C, 6.49 at 145 C". Eq. 22 gives **6.978 (→ 6.98), 6.467 (→ 6.47), 6.486 (→ 6.49)**.
  The Table 2 values it quotes (6.99 / 6.49 / 6.49 at 25 / 100 / 150 °C) are correct.
- The constant-ΔH fallback (`PKA_25C["h2s_1"] = 7.05`, ΔH 22.0 kJ/mol) gives [Z] 6.28 at 100 °C,
  **5.94 at 145 °C** and 5.62 at 200 °C, against 6.47 / **6.49** / 6.74 from Eq. 22: **0.54 pK
  units too acidic at 145 °C**. Its 25 °C anchor (7.05) sits 0.06 above Table 2's 6.99 (0.09 above
  the 100-bar 6.96). Its ΔH (22.0) is close to Table 3's 22.7 kJ/mol. The failure is the missing
  ΔCp (−265), not the ΔH.
- The function returns the I = 0 molal pK at s.v.p. Ionic strength is applied separately
  (`pka_at` subtracts `davies_term`), consistent with the paper's I < 0.03 conditions.

---

## 6. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Can anchor:**
- **pK1 of H2S from 25 to 350 °C** (molal, I → 0, s.v.p.), measured, ±0.013–0.023 to 200 °C.
  This covers the whole Maillard window (100–180 °C) by interpolation, which the authors permit.
  It replaces any van 't Hoff extrapolation.
- **ΔH and ΔCp of ionisation** at 25 °C (22.7 kJ/mol, −265 J mol⁻¹ K⁻¹; the latter imposed from
  Barbero 1982) and their temperature course (Table 3, 100 bar).
- Pressure is immaterial for food (≤ 0.03 in log K up to 300 °C).

**Cannot anchor:**
- **Second ionisation (HS⁻ → S²⁻):** not measured.
- **High ionic strength / food matrices:** I < 0.03 only. The I-correction at 0.5–1.5 M (Davies) is
  the model's own assumption.
- **Below 25 °C:** Eq. 22 is extrapolative there. Irrelevant for heating, but relevant if a cooled-
  pH label scale is ever computed at < 25 °C.
- **Volatility / Henry's law:** this paper uses Suleimenov & Krupp 1994's Henry constants for Table 4
  but does not measure them. Not transcribed here.
- **Kinetics:** none. This is an equilibrium constant only.

**Caveats.** A single laboratory and method. The authors note agreement with Barbero 1982 at low to
moderate T and divergence from it (~0.3 log units) above 200 °C; Fig. 8 shows the older data. At
food temperatures the spread among sources is small next to the 0.5-unit constant-ΔH error. The
[Z] numbers were computed in the `maillard` env from Eq. 22 and the model function's source. The
model was not run.

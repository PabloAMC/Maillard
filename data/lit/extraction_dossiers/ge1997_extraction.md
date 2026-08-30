# Ge & Lee 1997 (J. Agric. Food Chem. 45, 1619–1623; JF960458D / S0021-8561(96)00458-X) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/ge1997.pdf` (5 pp.). Read method: both — born-digital text layer
(clean) plus 200–220 dpi page rasters of PDF pp. 1, 3, 4 and 5 to verify the abstract, Table 1,
Tables 2 and 3, and Figure 6 respectively. Every table value below is raster-confirmed.

---

## 0. TWO PROMINENT WARNINGS, PLUS PAPER IDENTITY

### WARNING (a) — THE METHOD OF THIS PAPER IS DISPUTED IN PRINT BY VAN BOEKEL

`docs/reference/FIT_HOLDOUT_DECLARATION.md` line 218 carries the row:

> `| Ge & Lee 1997 reverse-rate Ea | **quarantined until tables retrieved** | method disputed in print by van Boekel; abstract has a unit defect (kcal/J confusion) |`

**I searched the repository for the content of that critique and it is NOT on disk.** Greps over
`data/` and `docs/` for `Ge 1997`, `ge1997`, `Ge and Lee`, `Ge & Lee`, `Ge, S`, `1619` and
`van Boekel` return exactly one hit — the declaration line quoted above — plus unrelated matches.
There is no note, no quotation, no citation, and no page reference anywhere in the repository that
states **what** van Boekel objected to, **where** he said it, or **which** part of the method it
touches.

**Therefore, stated as a flagged caveat and nothing more:** *the method of Ge & Lee 1997 is
disputed in print by van Boekel; the specific content of that critique is not recorded anywhere in
this repository and is NOT reproduced here, because inventing it would be fabrication.* Any wave
that wants to act on the dispute must first retrieve and cite the actual van Boekel text.

What this extraction *can* offer instead is a set of **independently verifiable, self-contained
defects in the paper's own arithmetic and internal consistency**, documented in §6 and §7 below,
which stand on their own regardless of what van Boekel wrote.

### WARNING (b) — THE ABSTRACT CARRIES A UNIT DEFECT, AND **SO DOES TABLE 3**; NEITHER "kcal/mol" NOR "J/mol" IS SELF-CONSISTENT

Fully worked in §7. Summary of the resolution, up front:

- The **abstract** and **Table 3** print the SAME numbers with the SAME unit label, `kcal/mol`. There
  is **no abstract-vs-table discrepancy** — the defect is present identically in both places.
- Read as printed (**kcal/mol**): 6.52 × 10⁴ kcal/mol = 2.73 × 10⁸ J/mol ≈ 273 MJ/mol. Physically
  absurd (a C–C bond is ~350 kJ/mol). **REJECT.**
- Read as **J/mol**: 6.52 × 10⁴ J/mol = 65.2 kJ/mol. **Also rejected** — the paper's own Table 2
  k₁ values give an Arrhenius slope of 118.4 kJ/mol, not 65.2.
- Read as **cal/mol**: 6.52 × 10⁴ cal/mol = 65.2 kcal/mol = 273 kJ/mol. **Also rejected** as a
  correct Ea — but this *is* the arithmetic route the authors took, and it is off by exactly a factor
  of 2.303 from the correct value.
- **The correct activation energies, recomputed from the authors' own Table 1 and Table 2 k-values
  (my regression reproduces their printed correlation coefficients to 4 decimal places), are
  Ea₁ = 28.3, Ea₋₁ = 32.4, Ea₂ = 34.8 kcal mol⁻¹ (= 118.4, 135.6, 145.5 kJ mol⁻¹).**
- **I have NOT silently picked one reading.** All three are shown with arithmetic in §7; none is
  self-consistent; the number is reproducible only as `slope(ln k vs −1/T) × 2.303R`, which
  double-counts the ln↔log₁₀ conversion factor.

---

## 0.1 PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY

| field | value |
|---|---|
| authors | Shi-Jun Ge and Tung-Ching Lee* (corresponding) |
| title | "Kinetic Significance of the Schiff Base Reversion in the Early-Stage Maillard Reaction of a Phenylalanine−Glucose Aqueous Model System" |
| venue | *Journal of Agricultural and Food Chemistry* |
| volume/pages/year | **45**, No. 5, 1619–1623, 1997 |
| DOI | not printed; manuscript ID **JF960458D**; CCC code S0021-8561(96)00458-X |
| received/revised/accepted | "Received for review June 24, 1996. Revised manuscript received February 7, 1997. Accepted February 18, 1997." Abstract published in *Advance ACS Abstracts*, April 1, 1997 |
| affiliation | Center of Advanced Food Technology, Department of Food Science, and Institute of Marine and Coastal Sciences, Cook College, Rutgers, The State University of New Jersey, New Brunswick, New Jersey 08903-0231. Publication D-10114-2-97 of the New Jersey Agricultural Experiment Station |
| PDF character | born-digital (Acrobat Distiller 2.1 for SunOS/Solaris, 1997), clean text layer, 5 pp. |

Expected identity was "Schiff base reversion kinetics". **Confirmed.**

## 1. ONE-PARAGRAPH VERDICT

**Carrying warnings (a) and (b) forward: the method of this paper is disputed in print by van Boekel
(content not on disk — see §0), and the reported activation energies are unit-defective in both the
abstract and Table 3 (see §7).** With those flags standing, the paper gives the model a
seven-temperature (40–97 °C) set of rate constants for a glucose + phenylalanine aqueous model system
at pH 7.0: the Amadori-rearrangement constant k₂ (h⁻¹) and equilibrium constant K in Table 1, and the
Schiff-base reversion constant k₋₁ (h⁻¹) plus the derived formation constant k₁ (mM⁻¹ h⁻¹) in
Table 2 — **21 rate constants and 7 equilibrium constants in total, all measured at a single pH.** It
also gives three activation energies (Table 3) which are numerically wrong (§7), and one
figure-only pH profile at 80 °C spanning pH 2–12 (Figure 6). It gives regression quality statistics
which are themselves mislabelled: the column headed `R²` is in fact **R**, verified to 4 decimal
places (§6.1). It gives **no error bars of any kind**, no per-point SD, and no n beyond "Triplicate
samples were taken for each concentration, and at least seven different concentrations were tested
for each temperature". The K column has **no printed units** and the k₁ column can only be
reconciled with K if K is read as M⁻¹ (§6.3), which the paper never states.

## 2. SYSTEM DEFINITION — verbatim

**Materials and Methods, p. 1620 (PDF page 2), verbatim:**
> "D-(+)-Glucose, sodium metabisulfite, ammonium hydroxide (29.1%), and sodium acetate were ACS
> analytical grade reagents from Sigma Chemical Co., St. Louis, MO. L-Phenylalanine was a SigmaUltra
> reagent. Amberlite CG120 (200-400 mesh), Dowex-50W (200-400 mesh), and Amberlite I-6766 (100-200
> mesh), were also purchased from Sigma. Other chemicals used in this study were all ACS analytical
> grade reagents.
> The model system was prepared by solubilizing equimolar ratios of glucose and phenylalanine, from
> 10 to 100 mmol, in a 0.05 M phosphate buffer, pH 7.0. A 15 mL scintillation vial was used as the
> reactor. Each vial was filled with 10 mL of reaction mixture and stored at −20 °C before use. The
> reaction was started when the temperature was raised to that of the water bath. The reaction was
> terminated by putting the vial into an ice-water bath and then transferring it to a −20 °C freezer
> for storage before compositional analysis.
> The phenylalanine Amadori compound, fructosylphenylalanine, was synthesized according to a
> procedure modified from that of Hashiba (1976) previously developed in our laboratory (Ge and Lee,
> 1996).
> The Maillard reaction mixture was analyzed using the HPLC method previously developed in our
> laboratory on a CarboPack PA-1 column, 4.6 × 250 mm (Dionex Corp., Sunnyvale, CA) (Ge and Lee,
> 1996). The HPLC system used for this study was a Dionex BioLC system (Dionex), which consisted of
> an advanced gradient pump, a UV detector, and a reagent storage module. The effective gradient
> program for the Maillard reaction mixture analysis was as follows: Elutant A (0.1 M ammonium
> hydroxide) was linearly decreased from 100% to 40%, whereas elutant B (0.1 M ammonium hydroxide-0.5
> M sodium acetate) was linearly increased from 0% to 40% in 16.7 min at the same time. A good
> separation of the reaction mixture using the HPLC method was obtained, and the chromatography is
> shown in Figure 1. The retention times for phenylalanine and fructosylphenylalanine were 11.73 and
> 13.5 min, respectively. The flow rate of the elutant was kept constant at 1 mL/min. Both
> phenylalanine and fructosylphenylalanine were detected using a UV detector at 260 nm. The data were
> acquired and recorded by a Waters chromatography workstation Maxima 820 (Waters Dynamic Solutions,
> Ventura, CA). A calibration curve was used for determination.
> Statistical analysis of the data was done using linear regression of the appropriate rate function
> listed in the above section to determine their rate constants."

| condition variable | value as printed |
|---|---|
| temperature | **40, 54, 68, 76, 81, 90, 97 °C** (7 levels, water bath). Figure 6 is at **80 °C** |
| buffer | **0.05 M phosphate buffer** for the kinetics. Figure 6 uses a buffer ladder: "0.01 N hydrochloric acid for pH 2.0, 0.05 M acetic buffer for pH 4.0, 0.05 M phosphate buffer for pH 6.0 and 8.0, 0.05 M boric buffer for pH 10, and 0.01 N sodium hydroxide for pH 12" |
| pH | **7.0** for all of Tables 1–3 and Figures 1–5. Figure 6 spans **2.0, 4.0, 6.0, 8.0, 10, 12** |
| hot-pH correction | **NOT DONE and not discussed.** pH 7.0 is a room-temperature buffer formulation; the reactions run at 40–97 °C. Phosphate buffer pH drifts with temperature and this is never addressed |
| concentrations | "equimolar ratios of glucose and phenylalanine, from **10 to 100 mmol**" (Table 1 footnote: "in the range of 10–100 mM"). Figures 1, 2, 6 use **0.05 M** each |
| vessel | 15 mL scintillation vial |
| volume | 10 mL per vial |
| headspace / atmosphere | 5 mL headspace implied; **atmosphere never stated**, no purge, no seal described |
| agitation | **none stated** |
| quench | ice-water bath, then −20 °C freezer |
| analytical method | HPLC, Dionex BioLC, CarboPack PA-1 4.6 × 250 mm, gradient of 0.1 M NH₄OH → 0.1 M NH₄OH + 0.5 M NaOAc over 16.7 min, 1 mL/min, UV 260 nm; t_R 11.73 min (Phe) and 13.5 min (fructosylphenylalanine) |
| calibration / internal standard | "A calibration curve was used for determination." **No internal standard.** No recovery data |
| replication (n) | "**Triplicate** samples were taken for each concentration, and **at least seven different concentrations** were tested for each temperature" (Table 1 footnote) |
| error bars | **NONE anywhere in the paper.** No ±, no SD, no SE, no confidence interval on any k, K or Ea |

**Theoretical model, p. 1619–1620 (PDF pages 1–2), verbatim.** Scheme:
`A + S ⇌(k₁, k₋₁) [AS]sb* →(k₂) AR`  (eq 1), "where A is amino acid, S is sugar, [AS]sb* is the
Schiff base complex intermediate, and AR is the Amadori compound."

Assumptions, verbatim: "We assume that at the early stages of the Maillard reaction (1) the
decomposition of amino acid through the Strecker degradation and of sugar through caramelization
will be negligible and (2) the conversion of the Amadori compound to the advanced Maillard reaction
products and its reversion to its parent compound or compounds will also be negligible at a low
concentration."

Rate equations as printed:
```
−d[A]/dt = k1[A][S] − k−1[AS]                              (2)
 d[AS]/dt = k1[A][S] − (k−1 + k2)[AS]                      (3)
 d[AR]/dt = k2[AS]                                         (4)
        K = [AS]/([A][S])                                  (5)
     [AS] = [A]0 − [A] − [AR]                              (6)
        K = ([A]0 − [A] − [AR]) / ([A][A])                 (7)
 d[AR]/dt = k2 K [A]²                                      (8)
−d[A]/dt = K k−1 [A]² − k−1([A]0 − [A] − [AR])             (9)
```
Determination procedure, verbatim: "For a given temperature and pH, the equilibrium constant K value
should not vary with the change of the concentrations of reactants. The K value can be calculated
from eq 7, and it should be equal to the average K value of different initial reactant
concentrations. The k₂ value can be calculated from eq 8 with the plot of d[AR]/dt against [A]². The
slope of the line equals the value of k₂K. The k₋₁ value can be calculated from eq 9 with the plot of
d[A]/dt against K[A]² − ([A]₀ − [A] − [AR]) of the different initial phenylalanine and glucose
concentrations at the early stage of the reaction. **The k₁ value can be calculated from the K and
k₋₁ values.**"

That last sentence fixes the provenance: **k₁ is derived, not measured.**

## 3. TABLE 1 — Amadori rearrangement k₂ and equilibrium constant K

**Anchor: Table 1, p. 1621 (PDF page 3).** Title as printed:
"Table 1. Linear Regression of the Experimental Data with Equation 8ᵃ"

Column headers as printed: `temp (°C)` | `k₂ (1/h)` | `R²ᵇ` | `K`
**Units: `1/h` for k₂; the K column carries NO units at all.**

Footnotes as printed:
> "ᵃ Conditions: concentrations of phenylalanine and glucose were equal and in the range of 10−100 mM
> in 0.05 M phosphate buffer solution, pH 7.0. Triplicate samples were taken for each concentration,
> and at least seven different concentrations were tested for each temperature.
> ᵇ The critical value for Spearman rank correlation coefficient is rs = 0.929 for α = 99%."

| temp (°C) | k₂ (1/h) | R²ᵇ | K | prov. |
|---|---|---|---|---|
| 40 | 1.27 × 10⁻³ | 0.9320 | 1.47 | [F] |
| 54 | 1.19 × 10⁻³ | 0.9605 | 1.09 | [F] |
| 68 | 6.74 × 10⁻² | 0.9772 | 0.71 | [F] |
| 76 | 2.04 × 10⁻¹ | 0.9979 | 0.53 | [F] |
| 81 | 2.29 × 10⁻¹ | 0.9415 | 0.96 | [F] |
| 90 | 1.007 | 0.9755 | 0.85 | [F] |
| 97 | 3.97 | 0.9887 | 0.37 | [F] |

Raster-verified against PDF p. 3; no cell unreadable. All values are **[F]** — regression slopes of
eq 8 (k₂) and averages over concentrations of eq 7 (K), not direct measurements.

**Defects in this table, reported not fixed:**
- **k₂ is non-monotonic in temperature:** 1.27 × 10⁻³ at 40 °C is *larger* than 1.19 × 10⁻³ at 54 °C.
  A rate constant that decreases with a 14 °C rise contradicts the paper's own Arrhenius treatment.
- **K is non-monotonic:** 1.47, 1.09, 0.71, 0.53, **0.96**, 0.85, 0.37. The 81 °C value breaks the
  otherwise decreasing trend by nearly a factor of 2.
- Two of the seven `R²` values (0.9320 at 40 °C and 0.9415 at 81 °C) are **below** the critical value
  0.929… no — 0.9320 and 0.9415 are both above 0.929, but only just. The text nevertheless claims
  "the correlation coefficients were much greater than the critical value 0.929", which is not true
  of those two rows. Both statements quoted here rather than reconciled.

## 4. TABLE 2 — Schiff base reversion k₋₁ and derived formation k₁

**Anchor: Table 2, p. 1622 (PDF page 4).** Title as printed:
"Table 2. Linear Regression of the Experimental Data with Equation 9ᵃ"

Column headers as printed: `temp (°C)` | `k₋₁ (1/h)` | `R²ᵇ` | `k₁ (mM⁻¹ h⁻¹)`
**Units: `1/h` for k₋₁; `mM⁻¹ h⁻¹` for k₁.**

Footnotes as printed:
> "ᵃ The experiment conditions were the same as those of Table 1.
> ᵇ The critical value for Spearman rank correlation coefficient is rs = 0.714 for α = 90%."

| temp (°C) | k₋₁ (1/h) | R²ᵇ | k₁ (mM⁻¹ h⁻¹) | prov. |
|---|---|---|---|---|
| 40 | 1.33 × 10⁻³ | 0.8942 | 1.97 × 10⁻⁶ | k₋₁ [F]; k₁ [F]/[Z] derived |
| 54 | 2.40 × 10⁻² | 0.6046 | 2.63 × 10⁻⁵ | as above |
| 68 | 1.29 × 10⁻¹ | 0.7398 | 1.54 × 10⁻⁴ | as above |
| 76 | 3.26 × 10⁻¹ | 0.7172 | 1.73 × 10⁻⁴ | as above |
| 81 | 4.35 × 10⁻¹ | 0.7417 | 4.90 × 10⁻⁴ | as above |
| 90 | 1.26 | 0.7150 | 1.07 × 10⁻³ | as above |
| 97 | 9.42 | 0.9484 | 3.54 × 10⁻³ | as above |

Raster-verified against PDF p. 4; no cell unreadable.

**Defects in this table, reported not fixed:**
- **The fit quality is poor.** Five of seven `R²` values are between 0.6046 and 0.7417. The 54 °C
  value (**0.6046**) is *below* the stated critical value of 0.714, yet the text asserts (p. 1621)
  "the data also fit well to eq 9 because all of the correlation coefficient values were >0.714 (α =
  90%)." **The table contradicts the text.** Both quoted.
- **k₁ is derived from K and k₋₁, not measured** (Methods, §2). It is therefore not an independent
  datum.
- **The `R²` column is actually R, not R²** — see §6.1.

## 5. TABLE 3 — activation energies (DEFECTIVE — see §7)

**Anchor: Table 3, p. 1622 (PDF page 4).** Title as printed:
"Table 3. Activation Energy (Ea) of the Maillard Reactions"

Header structure as printed: a spanning header `reaction step` over three columns `k₁` | `k₋₁` | `k₂`;
row labels `Eₐ (kcal/mol)` and `R²ᵃ`.
Footnote as printed: "ᵃ The critical value for Spearman rank correlation coefficient is rs = 0.929
for α = 99%."

| | k₁ | k₋₁ | k₂ | prov. |
|---|---|---|---|---|
| **Eₐ (kcal/mol)** | 6.52 × 10⁴ | 7.49 × 10⁴ | 8.01 × 10⁴ | [F] — **unit-defective, see §7** |
| **R²ᵃ** | 0.9908 | 0.9881 | 0.9685 | [F] — **mislabelled, is R not R², see §6.1** |

Raster-verified against PDF p. 4. **The unit label `kcal/mol` in Table 3 is identical to the abstract's.**

Supporting figures: **Figure 4** ("Effect of the reaction temperature on the kinetic constants k₁,
k₋₁, and k₂", p. 1622 / PDF p. 4) and **Figure 5** ("Linear regression plots of Ln(k) with −1/T",
p. 1622 / PDF p. 4). **Figure 5's y-axis is labelled `Ln(k)` and its x-axis `−1/T (10⁻³)`, spanning
−3.3 to −2.6, with Ln(k) from −14 to +4.** This is the load-bearing evidence for §7: **the regression
was on a NATURAL log**, so the Arrhenius slope is Eₐ/R, not Eₐ/(2.303R).

## 6. INTERNAL-CONSISTENCY CHECKS — three findings

### 6.1 The column headed `R²` is R, verified to 4 decimal places.

I regressed `ln k` on `(−1/T)` (T in K) using the paper's own Table 1 and Table 2 values, exactly as
Figure 5 depicts. Pearson correlation and its square:

| step | my r (Pearson) | my r² | printed as `R²` | verdict |
|---|---|---|---|---|
| k₁ | **0.99078** | 0.98164 | **0.9908** | printed value = **r**, not r² |
| k₋₁ | **0.98817** | 0.97648 | **0.9881** | printed value = **r** |
| k₂ | **0.96861** | 0.93820 | **0.9685** | printed value = **r** |

All three match to 4 decimal places on **r**, and none matches r². **The Table 3 `R²` header is
wrong; the values are correlation coefficients.** This also proves my regression reproduces the
authors' regression exactly, which is what licenses the Eₐ arithmetic in §7.

### 6.2 Table 2's `R²` values contradict the text (§4). The 54 °C value 0.6046 is below the stated 0.714 threshold.

### 6.3 The K column has no units, and its implied unit is M⁻¹ — but Table 1 and Table 2 disagree at two temperatures.

The Methods say "The k₁ value can be calculated from the K and k₋₁ values", i.e. k₁ = K · k₋₁.
Testing that with the printed columns (k₁ in mM⁻¹ h⁻¹, k₋₁ in h⁻¹):

| T (°C) | K (Table 1) | k₋₁ (Table 2) | K · k₋₁ | k₁ printed (mM⁻¹ h⁻¹) | (K·k₋₁)/k₁ | K implied by k₁ (= 1000·k₁/k₋₁) |
|---|---|---|---|---|---|---|
| 40 | 1.47 | 1.33 × 10⁻³ | 1.955 × 10⁻³ | 1.97 × 10⁻⁶ | 992 | 1.48 ✓ |
| 54 | 1.09 | 2.40 × 10⁻² | 2.616 × 10⁻² | 2.63 × 10⁻⁵ | 995 | 1.10 ✓ |
| 68 | 0.71 | 1.29 × 10⁻¹ | 9.159 × 10⁻² | 1.54 × 10⁻⁴ | **595** | **1.19 ✗** |
| 76 | 0.53 | 3.26 × 10⁻¹ | 1.728 × 10⁻¹ | 1.73 × 10⁻⁴ | 999 | 0.53 ✓ |
| 81 | 0.96 | 4.35 × 10⁻¹ | 4.176 × 10⁻¹ | 4.90 × 10⁻⁴ | **852** | **1.13 ✗** |
| 90 | 0.85 | 1.26 | 1.071 | 1.07 × 10⁻³ | 1001 | 0.85 ✓ |
| 97 | 0.37 | 9.42 | 3.485 | 3.54 × 10⁻³ | 984 | 0.376 ✓ |

All arithmetic in the table above is **[Z]**.

**Two conclusions:**
1. **The factor of ~1000 shows K is implicitly in M⁻¹ (l mol⁻¹), not mM⁻¹.** k₁ [mM⁻¹ h⁻¹] =
   K [M⁻¹] · k₋₁ [h⁻¹] / 1000. The paper never prints a unit for K; **M⁻¹ is inferred here, [Z].**
2. **At 68 °C and 81 °C, Table 1's K and Table 2's k₁ are mutually inconsistent** (implied K of 1.19
   and 1.13 vs printed 0.71 and 0.96). These are exactly the two rows that break the K monotonicity
   noted in §3. Reported, not fixed.

## 7. THE ACTIVATION-ENERGY UNIT DEFECT — full working

### 7.1 The two verbatim quotations

**ABSTRACT, p. 1619 (PDF page 1), raster-read, verbatim:**
> "The Eₐ₂ (8.01 × 10⁴ kcal/mol) was slightly greater than the Eₐ₁ (6.52 × 10⁴ kcal/mol) and the
> Eₐ₋₁ (7.49 × 10⁴ kcal/mol), revealing that the Amadori rearrangement was more sensitive to
> temperature changes and was more favorable at higher temperatures."

**TABLE 3, p. 1622 (PDF page 4), raster-read, verbatim (row label and cells):**
> `Eₐ (kcal/mol)` | `6.52 × 10⁴` | `7.49 × 10⁴` | `8.01 × 10⁴`

**Finding: the abstract and the table print IDENTICAL numbers with IDENTICAL units.** There is no
abstract-vs-table contradiction to adjudicate. The defect is in both places at once.

### 7.2 The three candidate readings, each tested against the paper's own k values

The reference calculation, using the paper's own Table 1 and Table 2 columns and the natural-log
Arrhenius form that Figure 5 depicts (`ln k = ln A + (Eₐ/R)·(−1/T)`, slope = Eₐ/R), with my
regression validated against their printed correlation coefficients in §6.1:

| step | slope of ln k vs (−1/T) | R = 8.314 J mol⁻¹ K⁻¹ → Eₐ | R = 1.987 cal mol⁻¹ K⁻¹ → Eₐ |
|---|---|---|---|
| k₁ | 14 240 K | **118 392 J/mol = 118.4 kJ/mol** | **28 298 cal/mol = 28.30 kcal/mol** |
| k₋₁ | 16 305 K | **135 555 J/mol = 135.6 kJ/mol** | **32 400 cal/mol = 32.40 kcal/mol** |
| k₂ | 17 498 K | **145 474 J/mol = 145.5 kJ/mol** | **34 771 cal/mol = 34.77 kcal/mol** |

All three slopes and both conversions are **[Z]**.

**Reading (a) — "kcal/mol" exactly as printed.**
6.52 × 10⁴ kcal/mol × 4184 J/kcal = **2.73 × 10⁸ J/mol = 273 MJ/mol**.
For scale, a C–C single bond is ~350 kJ/mol and the strongest bond in chemistry (CO) is ~1070
kJ/mol. 273 MJ/mol is ~250 000× a typical bond energy.
**REJECT — physically impossible, and 9650× the value the paper's own k data support.**

**Reading (b) — "J/mol".**
6.52 × 10⁴ J/mol = **65.2 kJ/mol**. The paper's own k₁ column gives **118.4 kJ/mol**.
Ratio = 118.4 / 65.2 = **1.816**. Not 1, not 2, not 2.303 — no clean factor.
**REJECT — not self-consistent with the paper's own data.**

**Reading (c) — "cal/mol".**
6.52 × 10⁴ cal/mol = **65.2 kcal/mol = 272.8 kJ/mol**. The paper's own k₁ column gives **28.30
kcal/mol**.
Ratio = 65.2 / 28.30 = **2.304**.
**REJECT as a correct Eₐ — but this ratio is exactly 2.303, and it identifies the arithmetic error.**

### 7.3 What the authors actually computed

`slope(ln k vs −1/T) × 2.303 × R(cal)`, i.e. `slope × 4.5766 cal mol⁻¹ K⁻¹`:

| step | slope × 4.5766 | printed in Table 3 and the abstract | agreement |
|---|---|---|---|
| k₁ | 14 240 × 4.5766 = **65 164** | 6.52 × 10⁴ | **0.06%** |
| k₋₁ | 16 305 × 4.5766 = **74 614** | 7.49 × 10⁴ | **0.4%** |
| k₂ | 17 498 × 4.5766 = **80 073** | 8.01 × 10⁴ | **0.03%** |

**All three reproduce to better than 0.5%.** The authors took the slope from a **natural-log** plot
(Figure 5, y-axis literally labelled `Ln(k)`) and multiplied it by **2.303 R**, the constant that
belongs to the **base-10** Arrhenius form `log k = log A − Eₐ/(2.303RT)`. **The 2.303 is spurious.**
The result is in **cal/mol**, and it was then labelled **kcal/mol**, adding a second, independent
factor-of-1000 error.

### 7.4 Resolution

**Two compounding defects, not one:**
1. a **log-base error** inflating every Eₐ by exactly 2.303, and
2. a **unit-prefix error** labelling cal/mol as kcal/mol (factor 1000).
Net: the printed numbers are **2303× too large** relative to a correctly computed Eₐ in kcal/mol.

**Neither "kcal/mol" nor "J/mol" rescues the printed numbers.** I have not picked one.

**The self-consistent values, recomputed from the paper's own tables ([Z], not printed by the
authors):**

| step | Eₐ (kcal mol⁻¹) [Z] | Eₐ (kJ mol⁻¹) [Z] |
|---|---|---|
| k₁ (Schiff formation) | **28.3** | **118.4** |
| k₋₁ (Schiff reversion) | **32.4** | **135.6** |
| k₂ (Amadori rearrangement) | **34.8** | **145.5** |

**One thing the defect does NOT damage:** the *ratios* are preserved. Printed
6.52 : 7.49 : 8.01 vs recomputed 28.30 : 32.40 : 34.77 → normalised to the same first term,
6.52 : 7.46 : 8.01. **The paper's qualitative conclusion — Eₐ₂ > Eₐ₋₁ > Eₐ₁, so the Amadori
rearrangement is the most temperature-sensitive step — survives intact.** Only the absolute values
and units are wrong.

## 8. THE PAPER'S OWN COMPARATIVE CLAIMS — and a dimensional problem in the headline

**Abstract, p. 1619, verbatim:** "The k₋₁ and k₂ were 10³ times greater than the k₁, indicating the
Schiff base formation, but not the Amadori rearrangement, was the rate-limiting step of the
reaction."

**Body, p. 1621, verbatim:** "The k₁ values were 10³ times smaller than both the values of k₋₁ and
k₂."

**[Z] check, and a defect:** k₁ is in **mM⁻¹ h⁻¹** and k₋₁, k₂ are in **h⁻¹**. **The comparison is
dimensionally invalid** — it compares a second-order rate constant to two first-order rate constants.
The numerical ratios also vary by a factor of 4 across the temperature range:

| T (°C) | k₋₁/k₁ | k₂/k₁ |
|---|---|---|
| 40 | 675 | 645 |
| 54 | 913 | 45 |
| 68 | 838 | 438 |
| 76 | 1884 | 1179 |
| 81 | 888 | 467 |
| 90 | 1178 | 941 |
| 97 | 2661 | 1121 |

All [Z]. The 54 °C k₂/k₁ ratio is **45**, not 10³. The "10³" claim is a loose average, not a
consistent finding, and it is dimensionally unsound as stated.

**Comparison to Higgins & Bunn, p. 1621, verbatim (a [C] datum):** "The results were similar to those
of the protein−sugar reaction system reported by Higgins and Bunn (1981) in the hemoglobin and
glucose reaction system, where **k₁ was 0.3 × 10⁻³ 1/mM·h and k₋₁ was 0.33 1/h** for the
glycosylation of the hemoglobin." — **[C]**, Higgins & Bunn, *J. Biol. Chem.* 1981, 256, 5204–5208.
(Note for cross-referencing: this is the same biomedical system as Weykamp & Penders 1982, whose
independently measured k₋₁ is 0.435 h⁻¹ — see `weykamp1982_extraction.md`. The two literature values
0.33 and 0.435 h⁻¹ agree to within 30%.)

## 9. FIGURE 6 — the pH profile (figure-only)

**Anchor: Figure 6, p. 1623 (PDF page 5).** Caption as printed, verbatim:
> "Figure 6. Effect of pH on the formation rate of the phenylalanine Amadori compound. Concentrations
> for phenylalanine and glucose were the same at 0.05 M. The reaction temperature was 80 °C. The pH
> buffer used was 0.01 N hydrochloric acid for pH 2.0, 0.05 M acetic buffer for pH 4.0, 0.05 M
> phosphate buffer for pH 6.0 and 8.0, 0.05 M boric buffer for pH 10, and 0.01 N sodium hydroxide for
> pH 12."

Axes: y = `dAR/dt (mM/h)`, 0 to 2 with ticks at 0, 1, 2; x = `pH`, 0 to 14 with ticks every 2.
Six data points connected by straight segments.

**Digitised from the 220-dpi raster of PDF page 5, tagged [fig] — approximate, ±0.05 mM/h:**

| pH | d[AR]/dt (mM/h) | provenance |
|---|---|---|
| 2.0 | ≈ 0.00 | **[fig]** |
| 4.0 | ≈ 0.00 | **[fig]** |
| 6.0 | ≈ 0.00 | **[fig]** |
| 8.0 | ≈ 0.05 | **[fig]** |
| 10 | ≈ 0.25 | **[fig]** |
| 12 | ≈ 1.67 | **[fig]** |

**These are the ONLY pH-dependent kinetic data in the paper, and they exist only as a figure.** No
table, no numbers in text. Note the pH ladder uses **five different buffer systems** (HCl, acetate,
phosphate, borate, NaOH) at two different formal concentrations, so ionic strength and specific-ion
catalysis are confounded with pH throughout this figure.

Text interpretation, verbatim, p. 1622–1623: "The rate of Amadori compound formation, as shown in
Figure 6, increased significantly with the increase of pH. This is due to the nucleophilic nature of
the Schiff base formation reaction. The active group of amino acid is amino group but not amine ion.
Under the alkaline pH conditions, the proton is released from the amine ion, increasing the effective
concentration of amino acid, which takes part in the reaction of the Schiff base formation."

## 10. OTHER FIGURES

| figure | anchor | content | numbers |
|---|---|---|---|
| Figure 1 (a, b) | p. 1621 / PDF p. 3 | "Typical HPLC chromatography of the Maillard reaction mixture of phenylalanine−glucose: (a) after development of the yellow-brown color; (b) before development of the yellow-brown color. Concentrations of phenylalanine and glucose were 0.05 M in 0.05 M phosphate buffer, pH 7.0. The samples were taken from the reaction mixture at the temperature of **76 °C for 24 h (a) and 8 h (b)**." Peak legend in (a): 1 phenylalanine; 2 fructosyl-phenylalanine; 3 brown color compound; 4 undetermined compound formed before peak 3 | t_R values only |
| Figure 2 | p. 1621 / PDF p. 3 | "Linear regression plots of d[AR]/dt with K[A]². Concentrations of the phenylalanine and glucose were 0.05 M in 0.05 M phosphate buffer, pH 7.0. The reaction temperatures ranged from 40 to 97 °C." Axes: d[AR]/dt (mM/h) −1 to 7; K[A]² 0 to 8. **Caption conflicts with Table 1's footnote**, which says the concentrations ranged 10–100 mM; the caption says 0.05 M. Both quoted | slopes = k₂K, already tabulated |
| Figure 3 | p. 1621 / PDF p. 3 | "Linear regression plots of −d[A]/dt with K[A]² − [A]₀ − [A] − [AR]. The reaction conditions were the same as those of Figure 2." Axes: −d[A]/dt (mM/h) 0 to 7; x 0 to 6 | slopes = k₋₁, already tabulated |
| Figure 4 | p. 1622 / PDF p. 4 | "Effect of the reaction temperature on the kinetic constants k₁, k₋₁, and k₂." Axes: Kinetical Constants 0 to 10; Reaction Temperature 30 to 100 °C. **Plots three quantities with two different unit systems on one unlabelled y-axis** | none beyond Tables 1–2 |
| Figure 5 | p. 1622 / PDF p. 4 | "Linear regression plots of Ln(k) with −1/T." Axes: **Ln(k)** −14 to 4; **−1/T (10⁻³)** −3.3 to −2.6. **This is the evidence that the regression was natural-log (§7.3)** | slopes reproduced in §7.2 |

## 11. DEFECTS AND CAVEATS TO CARRY FORWARD

1. **METHOD DISPUTED IN PRINT BY VAN BOEKEL. The content of that critique is not on disk in this
   repository and is not reproduced here.** (Warning (a), §0.)
2. **Eₐ unit defect, present identically in the abstract AND Table 3.** Neither kcal/mol nor J/mol is
   self-consistent; the numbers are cal/mol inflated 2.303× by a log-base error, hence 2303× too
   large as kcal/mol. Correct values: 28.3 / 32.4 / 34.8 kcal mol⁻¹ (118.4 / 135.6 / 145.5 kJ mol⁻¹).
   (Warning (b), §7.)
3. **The `R²` columns in Tables 1, 2 and 3 are R, not R²** — verified to 4 dp (§6.1).
4. **Table 2's text claim ("all >0.714") is falsified by its own 54 °C row (0.6046).**
5. **Table 1's text claim ("much greater than 0.929") is strained by its 40 °C (0.9320) and 81 °C
   (0.9415) rows.**
6. **k₂ is non-monotonic in T** (40 °C > 54 °C) and **K is non-monotonic** (81 °C breaks the trend).
7. **Table 1's K and Table 2's k₁ are mutually inconsistent at 68 °C and 81 °C** (§6.3).
8. **K has no printed units;** M⁻¹ is inferred here from the ×1000 factor, [Z].
9. **k₁ is derived (k₁ = K·k₋₁), not measured.** Not an independent datum.
10. **The "10³ times" headline compares a 2nd-order constant to two 1st-order constants** — a
    dimensionally invalid comparison — and the actual ratios span 45 to 2661 (§8).
11. **No error bars anywhere.** No SD, SE or CI on any k, K or Eₐ.
12. **No hot-pH correction.** pH 7.0 is a room-temperature formulation used at 40–97 °C.
13. **Atmosphere never stated; no agitation stated.**
14. **Figure 2's caption (0.05 M) contradicts Table 1's footnote (10–100 mM).**
15. **Figure 6's pH ladder confounds pH with five different buffer systems and two ionic strengths.**
16. The two "negligible" assumptions (no Strecker, no caramelisation, no Amadori degradation) are
    asserted and supported only by a chromatogram, never quantified.

## NEW-PARAMETER TABLE (consolidated)

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| k₂ (Amadori rearrangement) | 1.27 × 10⁻³ | 1/h | 40 °C, Phe+Glc equimolar 10–100 mM, 0.05 M phosphate, pH 7.0 | Table 1, p. 1621 / PDF p. 3 | [F] |
| k₂ | 1.19 × 10⁻³ | 1/h | 54 °C, as above | Table 1, p. 1621 / PDF p. 3 | [F] |
| k₂ | 6.74 × 10⁻² | 1/h | 68 °C, as above | Table 1, p. 1621 / PDF p. 3 | [F] |
| k₂ | 2.04 × 10⁻¹ | 1/h | 76 °C, as above | Table 1, p. 1621 / PDF p. 3 | [F] |
| k₂ | 2.29 × 10⁻¹ | 1/h | 81 °C, as above | Table 1, p. 1621 / PDF p. 3 | [F] |
| k₂ | 1.007 | 1/h | 90 °C, as above | Table 1, p. 1621 / PDF p. 3 | [F] |
| k₂ | 3.97 | 1/h | 97 °C, as above | Table 1, p. 1621 / PDF p. 3 | [F] |
| K (Schiff equilibrium) | 1.47 / 1.09 / 0.71 / 0.53 / 0.96 / 0.85 / 0.37 | **no units printed** (M⁻¹ inferred, §6.3) | 40 / 54 / 68 / 76 / 81 / 90 / 97 °C, pH 7.0 | Table 1, p. 1621 / PDF p. 3 | [F]; unit is **[Z]** |
| k₋₁ (Schiff reversion) | 1.33 × 10⁻³ | 1/h | 40 °C, pH 7.0 | Table 2, p. 1622 / PDF p. 4 | [F] |
| k₋₁ | 2.40 × 10⁻² | 1/h | 54 °C | Table 2, p. 1622 / PDF p. 4 | [F] |
| k₋₁ | 1.29 × 10⁻¹ | 1/h | 68 °C | Table 2, p. 1622 / PDF p. 4 | [F] |
| k₋₁ | 3.26 × 10⁻¹ | 1/h | 76 °C | Table 2, p. 1622 / PDF p. 4 | [F] |
| k₋₁ | 4.35 × 10⁻¹ | 1/h | 81 °C | Table 2, p. 1622 / PDF p. 4 | [F] |
| k₋₁ | 1.26 | 1/h | 90 °C | Table 2, p. 1622 / PDF p. 4 | [F] |
| k₋₁ | 9.42 | 1/h | 97 °C | Table 2, p. 1622 / PDF p. 4 | [F] |
| k₁ (Schiff formation) | 1.97 × 10⁻⁶ | mM⁻¹ h⁻¹ | 40 °C, pH 7.0 | Table 2, p. 1622 / PDF p. 4 | **[F]/[Z] — derived as K·k₋₁/1000, NOT measured** |
| k₁ | 2.63 × 10⁻⁵ | mM⁻¹ h⁻¹ | 54 °C | Table 2, p. 1622 / PDF p. 4 | as above |
| k₁ | 1.54 × 10⁻⁴ | mM⁻¹ h⁻¹ | 68 °C | Table 2, p. 1622 / PDF p. 4 | as above; **inconsistent with Table 1's K, §6.3** |
| k₁ | 1.73 × 10⁻⁴ | mM⁻¹ h⁻¹ | 76 °C | Table 2, p. 1622 / PDF p. 4 | as above |
| k₁ | 4.90 × 10⁻⁴ | mM⁻¹ h⁻¹ | 81 °C | Table 2, p. 1622 / PDF p. 4 | as above; **inconsistent with Table 1's K, §6.3** |
| k₁ | 1.07 × 10⁻³ | mM⁻¹ h⁻¹ | 90 °C | Table 2, p. 1622 / PDF p. 4 | as above |
| k₁ | 3.54 × 10⁻³ | mM⁻¹ h⁻¹ | 97 °C | Table 2, p. 1622 / PDF p. 4 | as above |
| **Eₐ₁ as printed** | 6.52 × 10⁴ | **kcal/mol (DEFECTIVE)** | 40–97 °C, pH 7.0 | Table 3 p. 1622 / PDF p. 4; abstract p. 1619 / PDF p. 1 | [F] — **DO NOT USE** |
| **Eₐ₋₁ as printed** | 7.49 × 10⁴ | **kcal/mol (DEFECTIVE)** | as above | Table 3 p. 1622 / PDF p. 4; abstract p. 1619 / PDF p. 1 | [F] — **DO NOT USE** |
| **Eₐ₂ as printed** | 8.01 × 10⁴ | **kcal/mol (DEFECTIVE)** | as above | Table 3 p. 1622 / PDF p. 4; abstract p. 1619 / PDF p. 1 | [F] — **DO NOT USE** |
| Eₐ₁ recomputed | 28.3 | kcal mol⁻¹ | 40–97 °C, pH 7.0 | §7.2, from Table 2 | **[Z]** |
| Eₐ₁ recomputed | 118.4 | kJ mol⁻¹ | as above | §7.2 | **[Z]** |
| Eₐ₋₁ recomputed | 32.4 | kcal mol⁻¹ | as above | §7.2, from Table 2 | **[Z]** |
| Eₐ₋₁ recomputed | 135.6 | kJ mol⁻¹ | as above | §7.2 | **[Z]** |
| Eₐ₂ recomputed | 34.8 | kcal mol⁻¹ | as above | §7.2, from Table 1 | **[Z]** |
| Eₐ₂ recomputed | 145.5 | kJ mol⁻¹ | as above | §7.2 | **[Z]** |
| Arrhenius pre-exponential slopes (Eₐ/R) | 14 240 / 16 305 / 17 498 | K | k₁ / k₋₁ / k₂, 40–97 °C | §7.2 | **[Z]** |
| correlation coefficient R (mislabelled `R²`) | 0.9908 / 0.9881 / 0.9685 | — | k₁ / k₋₁ / k₂ Arrhenius fits | Table 3, p. 1622 / PDF p. 4 | [F], **header defective** |
| d[AR]/dt vs pH at 80 °C | 0.00 / 0.00 / 0.00 / 0.05 / 0.25 / 1.67 | mM/h | pH 2 / 4 / 6 / 8 / 10 / 12, 0.05 M Phe + 0.05 M Glc, five different buffers | Figure 6, p. 1623 / PDF p. 5 | **[fig]** — digitised from the 220-dpi raster of PDF p. 5, ±0.05 |
| k₁ (haemoglobin–glucose), cited | 0.3 × 10⁻³ | 1/mM·h | 37 °C, biomedical | text p. 1621 / PDF p. 3 | **[C]** — Higgins & Bunn 1981 |
| k₋₁ (haemoglobin–glucose), cited | 0.33 | 1/h | 37 °C, biomedical | text p. 1621 / PDF p. 3 | **[C]** — Higgins & Bunn 1981 |
| replication | triplicate per concentration; ≥7 concentrations per temperature | — | — | Table 1 footnote a, p. 1621 / PDF p. 3 | [M] |
| error bars | **NONE PRINTED ANYWHERE** | — | — | — | — |

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

**Standing status.** `docs/reference/FIT_HOLDOUT_DECLARATION.md` line 218 currently reads
*"Ge & Lee 1997 reverse-rate Ea — **quarantined until tables retrieved** — method disputed in print
by van Boekel; abstract has a unit defect (kcal/J confusion)"*. The tables are now retrieved and
transcribed (§3, §4, §5). **My recommendation is that the quarantine be PARTIALLY lifted and
PARTIALLY hardened**, as follows — and that the lifting be conditional on someone first retrieving
the van Boekel text, because warning (a) is unresolved.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Table 3 Eₐ values as printed (6.52 / 7.49 / 8.01 × 10⁴ "kcal/mol")** | **PERMANENTLY EXCLUDED — HARDEN THE QUARANTINE** | — | Demonstrably wrong by a factor of 2303 (§7). No unit reading rescues them. These numbers must never enter any parameter, bound, initialisation or prior. |
| **Recomputed Eₐ: k₋₁ = 32.4 kcal/mol (135.6 kJ/mol)** | **HOLD-OUT — do not fit** | temperature | This is the one thing the paper uniquely offers the trunk: a temperature dependence for the Schiff-base **reversion**, the parameter that `weykamp1982` gives at 37 °C with an explicit `rate_transfer: not_licensed without the Ea` guard. **But it is [Z] — my recomputation, not the authors' printed value** — and it rests on a 7-point fit whose underlying k₋₁ regressions have R as low as 0.6046. Using it to *lift* the Weykamp guard would be circular in the worst way: importing a number I computed to license a transfer the declaration explicitly blocked. **Score against it; never fit to it.** |
| **Recomputed Eₐ: k₂ = 34.8 kcal/mol (145.5 kJ/mol), k₁ = 28.3 kcal/mol (118.4 kJ/mol)** | **HOLD-OUT — do not fit** | temperature | Same reasoning. Note k₁'s Eₐ is not independent: k₁ = K·k₋₁, so Eₐ₁ = Eₐ₋₁ + ΔH_K by construction. Two of the three "activation energies" are one measurement. |
| **Table 1 k₂ column, the 5 rows 68–97 °C** | **HOLD-OUT** | temperature (high arm) | The Arrhenius-clean part of the Amadori series. Excludes the 40 °C and 54 °C rows, which are non-monotonic (§3) and which drag the fit. |
| **Table 1 k₂, the 40 °C and 54 °C rows** | **EXCLUDED — internally inconsistent** | — | k₂(40 °C) > k₂(54 °C). A rate constant cannot fall with rising temperature under this mechanism. Do not score against either. |
| **Table 2 k₋₁ column, the rows at 40, 76, 90, 97 °C** | **HOLD-OUT** | temperature | The four rows whose K/k₁ cross-check closes to within 1% (§6.3). Cleanest subset of the reversion series. |
| **Table 2 k₋₁, the 54 °C row (R = 0.6046)** | **EXCLUDED — fails the authors' own stated threshold** | — | The paper claims all values exceed 0.714; this one does not. Excluding it is the conservative direction. |
| **Table 2 k₋₁, the 68 °C and 81 °C rows** | **EXCLUDED — Table 1 / Table 2 inconsistency** | — | At these two temperatures the printed K and the printed k₁ cannot both be right (§6.3, implied K 1.19 and 1.13 vs printed 0.71 and 0.96). One of the two tables is wrong at these rows and the print does not say which. |
| **Table 1 K column** | **DO NOT USE AS A THERMODYNAMIC K** | — | No printed units, non-monotonic in T, and inconsistent with k₁ at two of seven temperatures. Its unit (M⁻¹) is my inference, [Z]. If the trunk needs a Schiff-base K, the Weykamp K₁ = 2.12 M⁻¹ at 37 °C is unit-resolved and internally consistent to 4 digits; this one is neither. |
| **Table 2 k₁ column** | **DO NOT FIT SEPARATELY — DERIVED** | — | k₁ = K·k₋₁/1000 exactly (verified at 5 of 7 temperatures, §6.3). Fitting k₁ and k₋₁ and K would triple-count. |
| **Figure 6 pH profile (6 points, 80 °C, pH 2–12)** | **HOLD-OUT, low weight, [fig]-flagged** | pH | The only pH-resolved data in the paper, and orthogonal to everything else in it (all tables are pH 7.0 only). But: the values are digitised by me from a raster, three of six points sit on zero, and the pH ladder confounds pH with five different buffer systems and two ionic strengths. Usable only as a **shape/sign** test (rate rises steeply above pH 8, negligible below pH 6), never as a magnitude target. |
| **The "k₋₁ and k₂ are 10³ × k₁" claim** | **NOT A DATUM** | — | Dimensionally invalid (2nd-order vs 1st-order) and the actual ratios span 45–2661 (§8). |
| **Higgins & Bunn k₁ = 0.3 × 10⁻³ 1/mM·h, k₋₁ = 0.33 1/h** | **[C] — not a datum from this paper** | — | Cited, not measured here. If the trunk wants it, cite Higgins & Bunn 1981 directly. Note it is the same haemoglobin–glucose system as `weykamp1982` (k₋₁ = 0.435 h⁻¹); do **not** treat 0.33 and 0.435 as two independent constraints without checking whether Weykamp and Higgins & Bunn are independent measurements. |

**Circularity risks flagged — this paper carries several serious ones.**
(i) **The Eₐ values I recommend as HOLD-OUT are [Z], computed by me from the paper's own tables.**
They are not printed anywhere. Any wave that uses them must record that provenance. In particular,
using my recomputed Eₐ₋₁ to lift the `rate_transfer: not_licensed` guard on the Weykamp k₋₁ prior
would be justifying a transfer with a number this repository produced, not one the literature
printed. **I recommend against that specific move.**
(ii) **k₁ = K·k₋₁ exactly**, so Table 2's k₁ column and Table 1's K column and Table 2's k₋₁ column
are two independent quantities dressed as three. Likewise Eₐ₁ is not independent of Eₐ₋₁.
(iii) **The van Boekel dispute is unresolved (warning (a)).** Until its content is on disk, any FIT
role for this paper is being assigned in ignorance of a published objection to the method that
produced every number in it. **My recommendation is that nothing from this paper be given a FIT role
at all** — hold-out only — until that critique is retrieved and read.
(iv) **The `R²` mislabelling (§6.1) means the paper's fit quality is systematically worse than it
looks.** A reader taking 0.9908 as R² would believe 99% of variance explained; the true figure is
98.2% for the best fit and **60.5%** for the worst row in Table 2.
(v) **No error bars anywhere.** Any prior built on these values needs a deliberately wide sigma; the
paper supplies no basis for a narrow one.

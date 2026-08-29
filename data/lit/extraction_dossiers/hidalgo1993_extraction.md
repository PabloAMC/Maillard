# Hidalgo & Zamora 1993 — COMPLETE TRANSCRIPTION
### *"Non-enzymatic Browning and Fluorescence Development in a (E)-4,5-Epoxy-(E)-2-heptenal/Lysine Model System"*
### ⭐ **THE CORPUS'S ONLY Eₐ FOR A CARBONYL + LYSINE REACTION MEASURED FROM 4 °C UPWARDS — and it measures COLOUR, not adduct.** Wave **K6b**, 2026-08-30.

**Full extraction of every number in `data/articles/hidalgo1993.pdf`.**
Read-only wave. No repo file outside `data/lit/extraction_dossiers/` created or modified.

**Read first:** `research_round4_nulls.md` §A.1 (which named this paper as the one lysine
result in the Instituto de la Grasa Eₐ set) and `zamora2013_extraction.md` (same lab, same
estimator, twenty years later).

**Provenance:** **[M]** measured · **[C]** cited · **[Q]** qualitative · **[D]** digitized · **[Z]** derived here · **[!]** flag.

---

## 0. IDENTITY — **MATCHES EXACTLY. No mis-file.**

| field | value as printed |
|---|---|
| Authors | **Francisco J. Hidalgo, Rosario Zamora** |
| Venue | ***Journal of Food Science* 58 (3), 1993, pp. 667–670** |
| DOI | **10.1111/j.1365-2621.1993.tb04352.x** ✔ **exactly the expected DOI** |
| Affiliation | **Instituto de la Grasa y sus Derivados, CSIC, Avda. Padre García Tejero 4, 41012 Sevilla, Spain** |
| Keywords as printed | *nonenzymatic browning, fluorescence, heptenal-lysine, model-system* |
| Dates | *"Ms received 7/27/92; revised 9/30/92; accepted 12/22/92"* |
| Funding | **Comisión Interministerial de Ciencia y Tecnología of Spain, Project ALI 91-0409**; Prof. E. Vioque and Mr. J. L. Navarro thanked |
| PDF character | **4 pages**, **scanned + OCR'd, and the OCR text layer is heavily corrupted** (§0a), **Table 1** (one table), Figures 1–5 |

### 0a. ⚠️ **[!] THE TEXT LAYER IS UNRELIABLE — EVERY NUMBER BELOW WAS RE-VERIFIED ON A 400 dpi RASTER**
This is a 1993 scan. `pdftotext` renders the abstract's Eₐ as *"66.5 and SOKJ/mol"* (S for 5,
O for 0) and Table 1's pH column as *"3 / 5 / 6 / 6.5 / 7 / 7.5 / 6 / (blank) / (blank) / 12"*
with a `kB` value of **"6.12"** at pH 6.5. **Wave K6b rasterized page 668 at 400 dpi and read
Table 1 visually.** The corrected table is §2. **The OCR value "6.12" is actually "0.12"** —
a 51× error that would have inverted the pH trend. **Register: any automated ingest of this
PDF's text layer is unsafe; use §2 of this dossier.**

---

## 1. ⭐ THE ONE-PARAGRAPH ANSWER

**Two activation energies — 66.5 kJ/mol (colour) and 50 kJ/mol (fluorescence) — for
(E)-4,5-epoxy-(E)-2-heptenal reacting with free L-lysine at pH 7 over 4 / 25 / 40 / 80 °C.**
They are the corpus's only Eₐ pair for a lipid carbonyl + **lysine** (as opposed to
phenylalanine, asparagine or a mercaptan), and the only Instituto de la Grasa Eₐ measured from
**refrigeration temperature** upwards rather than from 80–120 °C.

**⚠️ AND THE OBSERVABLE IS THE PROBLEM.** The measured quantities are **CIELAB colour
difference ΔE** and **fluorescence intensity FI**. **Neither is a concentration of anything.**
Both are integrated, downstream, multi-step optical readouts of a polymerizing browning
system. **These are Eₐ of BROWNING, not of adduct formation** — a distinction the corpus must
enforce, because `meynier2004_extraction.md` §6b established, by direct measurement, that
**fluorescent-pigment intensity is NOT a proxy for covalent binding extent** (hexanal produced
3.5–7× more fluorophore than t-2-hexenal while producing ~7× less lysine loss).
**⇒ `hidalgo1993`'s 66.5 / 50 kJ/mol may NOT be substituted into an aldehyde–protein covalent
sink term. They belong on the browning lane.**

**A third, quieter result that is arguably the paper's most transferable: Table 1's ten-point
pH series, showing kB rising 116× and kF rising 864× from pH 3 to pH 12 — the corpus's densest
pH curve for a carbonyl–amine reaction.**

---

## 2. ⭐⭐ TABLE 1 — THE pH SERIES — TRANSCRIBED FROM A 400 dpi RASTER **[M]**

*Caption verbatim: "**Table 1** — Zero-order rate constant for browning (k_B) and fluorescence
(k_F) in an epoxyheptenal/lysine system as a funtion of pH ᵃ" (sic, "funtion").
Footnote a: "**Values given in hr⁻¹.** Experimental conditions are described in methods."*
*Conditions: **25 °C**, 0.1 mmol epoxyaldehyde + 0.2 mmol lysine in 5 mL buffer (§3).*

| **pH** | **k_B (hr⁻¹)** | **k_F (hr⁻¹)** | buffer |
|---:|---:|---:|---|
| **3** | **0.018** | **0.059** | 0.3 M sodium citrate |
| **5** | **0.029** | **0.041** | citrate |
| **6** | **0.059** | **0.80** | 0.3 M sodium phosphate |
| **6.5** | **0.12** | **2.30** | phosphate |
| **7** | **0.33** | **3.86** | phosphate |
| **7.5** | **0.41** | **5.29** | phosphate |
| **8** | **0.36** | **7.40** | 0.3 M sodium borate |
| **9** | **1.05** | **13.72** | borate |
| **10** | **1.53** | **22.10** | borate |
| **12** | **2.08** | **51.0** | 0.3 M sodium phosphate |

**[Z] Derived structure of the pH curve:**
- **k_B: pH 3 → 12 spans 116×.** Monotone **except** the pH 7.5 → 8 step (0.41 → 0.36, a 12 % dip)
  — **exactly the step where the buffer changes from phosphate to borate**. ⚠️ **[!] A buffer
  artefact is the parsimonious reading; borate is not inert toward polyols/amines.**
- **k_F: pH 3 → 12 spans 864×**, and is **monotone at every step** including the buffer change.
- **k_F/k_B rises from 3.3 (pH 3) to 24.5 (pH 12)** — the two observables do not scale together
  across pH, which is itself evidence they are not measuring the same species (cf. §5b).
- **The largest single jump for both is pH 5 → 6 (k_B 2.0×, k_F 19.5×).** [Z] The ε-NH₂ of
  lysine has pKₐ ≈ 10.5 and the α-NH₂ ≈ 9.0, so a 19.5× jump between pH 5 and 6 is **too
  steep and too low** to be free-base amine availability alone. Some of it is the
  citrate→phosphate buffer change. **Do not fit a pKₐ to this series.**
- Authors' summary [Q]: *"The data in Table 1 also suggest a **basic catalysis** for the
  browning reaction produced in the studied model system."*

**⭐⭐ AND THE PAPER'S OWN WARNING THAT RATE ≠ EXTENT, verbatim [M]:**
> *"At the end of the period (**20 days**) the higher brown color was reached at **neutral** pHs.
> However, the **browning rate was higher at higher pHs**. This apparent contradiction was due
> to a **higher initial reaction rate at higher pHs, but at these pHs the color formation
> stopped earlier**."*
And for fluorescence: *"After 20 days, higher fluorescence was obtained at **pH around 7 and
10**. However, the fluorescence rate was again higher at higher pHs."*
**⇒ RATE AND 20-DAY EXTENT HAVE DIFFERENT pH OPTIMA. Register
`pH_optimum_rate: monotone_increasing` and `pH_optimum_extent_20d: ~7 (colour), ~7 and ~10
(fluorescence)` as SEPARATE fields. This is the same rate/extent decoupling that
`shepelev2024_extraction.md` §4d finds for temperature, in a different system, in a different
decade — and it is the general failure mode of reading a kinetic constant as a yield.**

---

## 3. SYSTEM AND METHOD — complete transcription **[M]**

| variable | value as printed |
|---|---|
| **Carbonyl** | **(E)-4,5-epoxy-(E)-2-heptenal, 0.1 mmol** (standard), prepared from **(E)-2,(E)-4-heptadienal** per **Swoboda & Peers (1978)** |
| **Amine** | **L-lysine, 0.2 mmol** (standard) — **free amino acid, NOT a protein** |
| **Solvent** | **5 mL of 0.3 M sodium phosphate, pH 7.0** (standard) |
| **Dispersion** | *"suspended … and **sonicated to a stable emulsion** with a **Braun Labsonic U** sonicator"* ⚠️ **the epoxyaldehyde is not soluble — the system is an EMULSION, not a solution [!]** |
| **Standard temperature** | **25 °C** |
| **Temperature series** | **4, 25, 40 and 80 °C** ⭐ *(four points — the Arrhenius fit rests on four)* |
| **pH series** | citrate pH 3–5 · phosphate pH 6–7.5 · borate pH 8–10 · phosphate pH 12 |
| **Ratio series** | **0.1 and 0.2 mmol** epoxyaldehyde × **0.1, 0.2 and 0.4 mmol** lysine ⇒ E/L = **1:4, 1:2, 1:1, 2:1** |
| **Why this system [Q]** | the epoxyaldehyde is *"the **major volatile compound formed by the copper and α-tocopherol induced oxidation of butterfat**"* (Swoboda & Peers 1978) and a product of **ω-3 pentaenoic fatty acid** oxidation; lysine because *"it is usually lost during the deterioration of foods produced by peroxidizing lipids"* |

### 3a. Colour measurement **[M]**
200 µL of reaction + **2.8 mL deionized water** (a **15× dilution**). **Weighted-ordinate
method** (Hunter 1973). Transmittance on a **Hewlett-Packard 8450-A UV/VIS**, recorded at
**10 nm intervals from 400 to 700 nm**, **1 cm glass cells**; X,Y,Z → CIELAB L\*a\*b\* (CIE 1978).

- **ΔE = √[(Δa\*)² + (Δb\*)² + (ΔL\*)²]** vs a non-incubated control (Eq. 1)
- **YI_a = [100(1.28X − 1.06Z)]/Y** (Hunter 1942) (Eq. 2)
- **YI_b = 142.86 b\*/L\*** (Francis & Clydesdale 1975) (Eq. 3)

### 3b. Fluorescence measurement **[M]**
**Perkin-Elmer LS-5**, **100 µL** of reaction + **2.9 mL deionized water** (a **30× dilution** —
⚠️ **[!] a different dilution from the colour measurement, so ΔE and FI are not measured on the
same solution**). **Slit width 5 nm.** Standardized with **quinine sulfate (0.1 µM in 0.1 N
H₂SO₄) to give FI = 50 at 450 nm with excitation at 350 nm**.
**⇒ λ_ex 350 nm / λ_em 450 nm.** *(Compare `meynier2004_extraction.md` §6c: alkanal–amine
pigments 350/410, alkenal pigments 360–370/450. **This alkenal-derived system's 350/450 is a
hybrid of the two and is a calibration setting rather than a measured maximum — do not ingest
it as a spectral assignment.**)*

---

## 4. ⭐ THE ACTIVATION ENERGIES **[M]**

### 4a. The zero-order rate law, verbatim
Browning followed **zero-order kinetics (r = 0.992, p < 0.001)**:
> **ΔE = ΔE₀ + kt** (Eq. 4)

Yellowness index likewise (**r > 0.994, p < 0.001**):
> **YI = YI₀ + kt** (Eq. 5)

Fluorescence likewise (**r = 0.993, p < 0.001**):
> **FI = FI₀ + kt** (Eq. 6)

**Values at pH 7.0, 25 °C [M]:**

| observable | **intercept** | **k (hr⁻¹)** |
|---|---:|---:|
| **ΔE** | **0.4** | **0.328** |
| **YI_a** | **1.8** | **0.554** |
| **YI_b** | **1.2** | **0.483** |
| **FI** | **5** | **4.2** |

⚠️ **[!] THE TEXT'S k_FI = 4.2 hr⁻¹ DISAGREES WITH TABLE 1's 3.86 hr⁻¹ AT THE SAME pH 7 AND
25 °C** — an **8.8 %** internal inconsistency between the paper's own two reports of the same
number. (The ΔE values agree: text 0.328 vs Table 1 0.33 ✔.) **Register `k_FI(pH7,25C) =
3.86 or 4.2 hr⁻¹, paper-internal disagreement ±9 %`.**

### 4b. **THE Eₐ VALUES — verbatim from the Abstract and from Results**
> Abstract: *"Activation energy, according to the Arrhenius equation, was **66.5 and 50 KJ/mol
> for color difference and fluorescence intensity**, respectively."*
> Results: *"The rate constant values at different temperatures were calculated from Eq. (4)
> and (6), and used in an **Arrhenius plot (Fig. 4)** for calculation of activation energy (Eₐ)
> of color and fluorescence development. The values obtained for Eₐ of ΔE and FI were
> **66.5 and 50 KJ/mol**, respectively."*

| observable | **Eₐ (kJ/mol)** | n temperatures | range |
|---|---:|---:|---|
| **colour difference ΔE** | **66.5** | **4** | **4 → 80 °C** |
| **fluorescence intensity FI** | **50** | **4** | **4 → 80 °C** |

**⚠️ [!] "50" IS PRINTED WITHOUT A DECIMAL WHILE "66.5" CARRIES ONE.** The fluorescence value
is quoted to two significant figures at most; treat it as **50 ± 5**, not 50.0.
**No R², no confidence interval and no standard error is given for either fit.** Figure 4's
x-axis is labelled **1000/T** and its y-axis extends to about **−6**; the four points per
series are visible but **this wave does not digitize them, because the printed Eₐ values are
the paper's own fit and no useful precision would be added.**

### 4c. **The authors' own comparison, and it is a [C] not an [M]**
> *"These values were **one half of Eₐ value obtained for color changes in milk due to heat**
> (Pagliarini et al., 1990; and references therein), suggesting that **browning may be more
> easily produced by the model system we used**."*
**⇒ [C] implied Eₐ(milk colour on heating) ≈ 100–133 kJ/mol, from Pagliarini, Vernile & Peri
1990, *J. Food Sci.* 55, 1766–1767. Not measured here. Registered as a lead, not a value.**

### 4d. **[Z] What Eₐ = 66.5 / 50 kJ/mol actually implies, for the record**
Acceleration from 25 °C:

| Eₐ | 25 → 80 °C | 25 → 140 °C | 25 → 180 °C |
|---:|---:|---:|---:|
| **50** | **24×** | **352×** | **1 585×** |
| **66.5** | **73×** | **1 990×** | **13 900×** |

**⇒ These are large activation energies — the browning lane genuinely does accelerate steeply
with temperature. That is precisely why the contrast with
`shepelev2024_extraction.md`'s 15–20 kJ/mol for direct aldehyde–protein adduction is
informative rather than contradictory: browning and adduction are different reactions with
different rate-determining steps, and this paper's own §5b says so.**

---

## 5. THE OTHER MEASURED DEPENDENCES

### 5a. Epoxyaldehyde/lysine ratio — Figure 5 **[M]**
> *"An increase in the concentration of **either** epoxyaldehyde **or** lysine increased both
> color and fluorescence."*

| **E/L ratio** | **k for ΔE (hr⁻¹)** | **k for FI (hr⁻¹)** |
|---|---:|---:|
| **1 : 4** | **0.57** | **5.52** |
| **1 : 2** | **0.326** | **3.86** |
| **1 : 1** | **0.179** | **2.06** |
| **2 : 1** | **0.636** | **2.62** |

**[Z] Read carefully: the ratio series is NON-MONOTONE and has a MINIMUM at 1:1.**
Both observables fall from 1:4 → 1:1 and then rise again at 2:1. *(The 1:2 row is the standard
condition — 0.1 mmol aldehyde, 0.2 mmol lysine — and its values 0.326 / 3.86 reproduce Table 1's
pH 7 row 0.33 / 3.86 exactly ✔, a good internal consistency check.)*
Authors' own analysis [M]: *"The k values correlated with the **lysine concentration at the
same epoxyaldehyde concentration** for ΔE (**r = 0.9990, p < 0.5**), but **not for FI
(r = 0.98)**."*
⚠️ **[!] "p < 0.5" is not a significance statement** (and a correlation on **three** points at
fixed aldehyde is not either). **Do not ingest the reaction orders from this figure.**

### 5b. ⭐ The correlation between colour and fluorescence — and why it does not license substitution
Measured correlations [M]: ΔE ↔ YI ↔ FI **r = 0.9993, p < 0.001** as a function of time;
**r > 0.991, p < 0.001 at pH 5–8**; **r > 0.93, p < 0.01 at pH 9–12**; **r > 0.90, p < 0.001**
across temperature; **r > 0.98** across E/L ratio.

Conclusion, verbatim: *"This system produced color and fluorescence, **which always correlated**
in assayed conditions … Thus **the same reactions responsible for brown macromolecular pigment
formation may be related to production of fluorescent products**."*

**⚠️ [!] BUT THE TWO OBSERVABLES HAVE DIFFERENT ACTIVATION ENERGIES (66.5 vs 50 kJ/mol,
§4b), a k_F/k_B ratio that varies 7.4× across pH (§2), and different 20-day pH optima (§2).
Correlated ≠ identical. Three independent internal measurements in this very paper say the
two observables are NOT the same reaction coordinate.** And `meynier2004_extraction.md` §6b
supplies the decisive external falsification: in a single experiment, fluorophore intensity
and residue loss moved in **opposite** directions by ~7× each. **Register
`browning_and_fluorescence: correlated_within_this_system_only, NOT_interchangeable`.**

### 5c. Time course — Figure 1 **[M][Q]**
*"The reaction between epoxyaldehyde and lysine produced **significant changes in color and
fluorescence after 3 hr** of reaction."* Figure 1 runs **0–8 h at pH 7, 25 °C**, points are
**averages of triplicate determinations** — ⭐ **the only replication statement in the paper,
and it applies only to Fig. 1.** Figures 2, 3 and 5 and Table 1 carry **no n and no error bars**.

---

## 6. WHAT THIS PAPER DOES AND DOES NOT SUPPLY THE REPO

| question | answer |
|---|---|
| An Eₐ for lipid-carbonyl + **lysine** | ✅ **66.5 (colour) / 50 (fluorescence) kJ/mol, 4–80 °C** |
| An Eₐ for **adduct concentration** | ❌ **the observables are optical, not chemical (§1, §5b)** |
| A dense pH curve for a carbonyl–amine reaction | ✅ **ten points, §2 — the corpus's densest** |
| Ambient-to-moderate temperature coverage | ✅ **4, 25, 40, 80 °C — the only Instituto de la Grasa set that starts at refrigeration** |
| A protein | ❌ free L-lysine only |
| Reversibility | ❌ not tested |
| Saturated n-alkanal behaviour | ❌ the carbonyl is an epoxy-alkenal, a highly activated tertiary lipid oxidation product |
| Statistical rigour | ⚠️ n stated only for Fig. 1; no Eₐ error bars; one *"p < 0.5"* |

**SI: none. Nothing further to retrieve. The only outstanding lead is [C] Pagliarini et al.
1990 (milk colour Eₐ), §4c — low priority.**

---

## 7. INTEGRITY FLAGS

| # | flag | severity |
|---|---|---|
| **1** | **The OCR text layer corrupts Table 1 (kB at pH 6.5 reads "6.12", true value "0.12" — a 51× error) and the abstract ("SOKJ/mol" for "50 kJ/mol"). Use §2 of this dossier, never the text layer.** | ⭐⭐ **CRITICAL** |
| **2** | **Eₐ = 66.5 / 50 are BROWNING and FLUORESCENCE activation energies, not adduct-formation activation energies. `meynier2004_extraction.md` §6b falsifies fluorescence as a covalent-extent proxy by direct measurement.** | ⭐ HIGH |
| **3** | Colour and fluorescence are measured at **different dilutions** (15× vs 30×) on **different aliquots** | MEDIUM |
| 4 | k_FI at pH 7/25 °C is 3.86 (Table 1) or 4.2 (text) — a 9 % paper-internal disagreement | MEDIUM |
| 5 | The pH 7.5 → 8 dip in k_B coincides exactly with the phosphate → borate buffer change | MEDIUM |
| 6 | The reactant is a **sonicated emulsion**, not a solution; rates may carry an interfacial-area term | MEDIUM |
| 7 | *"r = 0.9990, p < 0.5"* on three points — not a significance statement | LOW–MEDIUM |
| 8 | Eₐ rests on **4** temperatures with no R², no SE, no CI; "50" is given to 2 s.f. | MEDIUM |
| 9 | n stated only for Fig. 1 (triplicate); absent everywhere else including Table 1 | MEDIUM |
| 10 | 350/450 nm is a **quinine-sulfate calibration setting**, not a measured fluorophore maximum | LOW |

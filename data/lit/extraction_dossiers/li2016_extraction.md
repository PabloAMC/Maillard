# Li, Jongberg, Andersen, Davies & Lund 2016 (10.1016/j.freeradbiomed.2016.05.019) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/li2016.pdf` (38 pp.). Read method: both — born-digital text layer
(`pdftotext -layout`) plus 150-dpi page rasters of PDF pp. 15, 21, 23, 38 to verify Tables 1, 2, 3
and the end of the file.

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY

| field | value |
|---|---|
| authors | Yuting Li (a,b,c), Sisse Jongberg (b), Mogens L. Andersen (b), Michael J. Davies (d), Marianne N. Lund (b,d, corresponding) |
| title | "Quinone-induced protein modifications: Kinetic preference for reaction of 1,2-benzoquinones with thiol groups in proteins" |
| venue | *Free Radical Biology and Medicine* — **Author's Accepted Manuscript** (unedited, pre-typesetting) |
| volume/pages/year | not yet assigned in this file; PII S0891-5849(16)30255-6, Reference FRB12882, 2016 |
| DOI | 10.1016/j.freeradbiomed.2016.05.019 |
| received/revised/accepted | Received 10 March 2016; Revised 11 May 2016; Accepted 18 May 2016 |
| affiliation | (a) School of Food Science and Technology, South China Univ. of Technology, Guangzhou; (b) Dept. of Food Science, Univ. of Copenhagen, Frederiksberg; (c) Guangdong Province Key Lab for Green Processing of Natural Products and Product Safety, SCUT; (d) Dept. of Biomedical Sciences, Univ. of Copenhagen |
| PDF character | born-digital, clean text layer, "ACCEPTED MANUSCRIPT" banner + "Accepted manuscript" watermark on every page |

Expected identity was "1,2-benzoquinone reacting with protein thiols; second-order rate constants
k2". **Confirmed.** The 1,2-benzoquinone studied is specifically **4-methyl-1,2-benzoquinone (4MBQ)**,
generated electrochemically from 4-methylcatechol (4MC).

**PAGINATION NOTE for all anchors below.** PDF p. 1 is the Elsevier "Author's Accepted Manuscript"
cover sheet. Manuscript page numbers are printed as footers, so **manuscript page N = PDF page N+1**.
Every anchor gives both.

**SUPPLEMENTARY INFORMATION IS NOT IN THIS PDF.** The paper cites "Tables 2–7, supplementary
material" (thiols), "Table 8, supplementary material" (Gly), "Table 9, supplementary material"
(Lys), and "Figure 1, supplementary materials" (cyclic voltammogram). PDF p. 38 (= ms p. 37) is the
TOC graphic + HIGHLIGHTS and is the **last page of the file** — verified by raster. None of the SI
tables are present. Every SI-derived number below is therefore recovered only from the way the
**main text quotes it**, and is marked as such.

## 1. ONE-PARAGRAPH VERDICT

This paper gives the model a complete, directly usable set of **apparent second-order rate constants
k2 (M⁻¹ s⁻¹) at 25 °C in 0.2 M phosphate buffer** for the reaction of a model o-quinone (4MBQ) with
three proteins (BSA, HSA, α-lactalbumin), one thiol-blocked protein (NEM-BSA), three low-molecular-
mass thiols (Nα-acetyl-L-Cys, L-Cys, GSH — as the "A" parameter of a saturating rate law), four
amines (Gly, L-Lys, Nε-acetyl-L-Lys, Nα-acetyl-L-Lys) and one guanidine (Nα-acetyl-L-Arg, an upper
bound only). It also gives an explicit **pH dependence** for the amines (pH 7.0 vs 8.0, Table 3),
for BSA (pH 7.0 vs 6.5, Table 1), and a *stated* pH for each low-molecular-mass thiol A-value
(pH 6.5 or 6.0). It gives **equilibrium constants K and unimolecular decay rate constants k
(s⁻¹)** for the tetrahedral quinone–thiol intermediate (Table 2). It gives cited thiol pKa values
(Nα-acetyl-L-Cys 9.5, L-Cys 8.3, GSH 8.8, all [C] from ref. 45). **There is NO Ea and NO temperature
series: every kinetic measurement in this paper is at a single temperature, 25 °C.** There is also
no ionic-strength series and no error model beyond the ± values printed (which the paper never
defines as SD or SE, and for which n is never stated). The full pH-resolved kobs tables that would
let one reconstruct the A(pH) and B(pH) curves are in a supplement that is not in this file.

## 2. SYSTEM DEFINITION — verbatim

**Quinone generation.** "On each day of analysis, a solution of 2 mM 4MC was dissolved in phosphate
buffer solution (0.2 M, pH 4.5) and de-oxygenated by purging with nitrogen for 10 min." Cyclic
voltammetry: CV-50W analyser (BAS), glassy carbon working electrode (3 mm dia.), platinum coil
counter electrode (5 mm dia., 5 cm long), Ag/AgCl (KCl, c = 3 M) reference (Metrohm). "The bulk
electrolysis was performed at an initial potential of 460 mV vs. Ag/AgCl under nitrogen … and the
initial concentration of 4MC solution was 2 mM. The working electrode was a reticulated vitreous
carbon tube with a diameter of 3.5 cm and length of 4.5 cm … Electrolysis was terminated when the
amount of charge passed was equivalent to 2 electrons per molecule. The 4MBQ yield was estimated to
be ca. 80% (1.6 mM) by UV-Vis spectrophotometry … using an extinction coefficient ε395 of 1350 M⁻¹
cm⁻¹ [28]. The 4MBQ solution was used within 1 h of generation, with the extent of loss being ca.
16% over this time" (ms p. 5–6 / PDF pp. 6–7).

**Kinetic experiments (the load-bearing methods paragraph, ms p. 7 / PDF p. 8), verbatim:**
> "The kinetic measurements were carried out at 25 °C by using a SX20 stopped-flow spectrophotometer
> (Applied Photophysics, London, UK) in the UV/Vis range (200-750 nm). Reactants were dissolved into
> 0.2 M phosphate buffer at appropriate pH values with 4MBQ in one syringe and the proteins, peptides
> or free amino acids in the other syringe. Equal volumes of each reactant were injected into the
> optical cell and UV-visible spectra recorded over time. The proteins, peptides or free amino acids
> were kept in at least 5-fold excess of 4MBQ concentration to obtain pseudo first-order kinetics."

| condition variable | value as printed |
|---|---|
| temperature | **25 °C** (kinetics); 22 °C for the amine-adduct preparation and Figure 1/Figure 9 spectra; 25 °C for LC column and DTNB incubation |
| buffer | **0.2 M phosphate** throughout the kinetics (0.2 M phosphate pH 4.5 for 4MBQ generation) |
| pH | 4.5 (generation) → mixed in the stopped flow with higher-pH phosphate to give **4.5–7.0** for the thiol/amine kinetic series; specific values used: 5.0, 5.5, 6.0, 6.5, 7.0, 8.0 (see per-table anchors). "the kinetic characterisation of the reactions of 4MBQ with thiols and amines in the range of pH values from 4.5 to 7.0" (ms p. 9 / PDF p. 10) — note this range statement **contradicts Table 3, which reports pH 8.0 data**; both quoted here rather than reconciled |
| hot-pH correction | not applicable (25 °C); no pH correction discussed |
| [4MBQ] | 0.04 mM in all stopped-flow runs; 0.2 mM for adduct-colour work; 1 mM for LC-MS/MS adduct prep |
| [protein] | HSA 0.2–0.5 mM (Fig. 3); BSA 0.28 mM and 0.5 mM (Figs. 4, 11); α-lactalbumin 0.28 mM; NEM-BSA 0.28 mM |
| [thiol] (LMM) | Nα-acetyl-L-Cys 0.5 mM (Fig. 5) and up to ~6 × 10⁻³ M in the concentration series (Fig. 6 x-axis); 10 mM for LC-MS/MS prep |
| [amine] | Nα-acetyl-L-Arg 50 mM; L-Lys 2.5 mM (Fig. 8); 0.5 mM (Fig. 11); 25 mM (Fig. 9 adduct spectra); 50 mM Gly / 100 mM Lys derivatives for LC-MS prep |
| vessel / volume | stopped-flow optical cell, equal volumes from two syringes; volumes not stated |
| headspace / atmosphere | N₂ purge during 4MBQ generation only; the stopped-flow runs are **not** stated to be deoxygenated. O₂ is explicitly required for the *coloured* amine-quinone product ("a fast re-coloration was observed in the presence of O₂") |
| agitation | none stated (stopped-flow); "stirred for 1 min" for adduct preparations |
| quench | none — continuous spectrophotometric monitoring |
| analytical method | UV-Vis stopped flow at 401 nm (4MBQ loss), 294 nm (thiol adduct), 499 nm (Nα-acetyl-L-Lys adduct), 508 nm (Nε-acetyl-L-Lys and L-Lys adducts), 293 nm; LC-ESI-MS/MS on ZIC-HILIC 2.1 mm × 15.0 cm, 3.5 µm, 0.1 mL/min, positive mode |
| thiol quantification | DTNB (Ellman), "based on the method of Jongberg et al. [23]": 600 µL protein + 2400 µL 0.1 M Tris pH 8.0 + 600 µL 10 mM DTNB, dark, 30 min, 25 °C, A₄₁₂ |
| calibration / internal standard | "Thiol concentrations were calculated based on a standard curve of L-Cys (1.67-83.3 µM) and given as mol thiol/mol protein"; ε395 = 1350 M⁻¹ cm⁻¹ for 4MBQ; no internal standard for LC-MS |
| NEM blocking | "The free thiol group in BSA (0.8 mM) was blocked by reaction with 80 mM NEM in 0.2 M phosphate buffer (pH 7.0) for 90 min at 25 C. Excess NEM was removed by eluting the protein through a PD-10 column (GE healthcare)." |
| replication (n) | **NEVER STATED.** No n anywhere in the paper |
| error bars | **NEVER DEFINED.** The paper prints "±" values in Tables 1, 2, 3 and in the text but never says whether they are SD, SE, or fit standard error. "Analysis of variance (ANOVA) was carried out by Microsoft Excel 2007. Significance was defined as a p value <0.05." is the entire statistics section |

**Rate law used for the LMM thiols (ms pp. 16–17 / PDF pp. 17–18), verbatim:**
> "Plots of kobs against the thiol concentrations did not give linear correlations. Similar non-linear
> correlations have been observed by Jameson et al. [30] for the reaction of dopaminoquione with Cys;
> these data were therefore analysed using equation 1, where [R-SH]T is the pseudo first-order
> concentration of thiol. … Equation 1 is based on a mechanism that involves an initial fast
> reversible formation of a tetrahedral intermediate between the quinone and the thiol, characterized
> by an equilibrium constant, K (see [30]). This tetrahedral intermediate decays by a first-order
> reaction with rate constant, k, to give the thiol-quinone adduct (Scheme 1)."

**Equations 1, 2 and 3 are images in this manuscript and carry no extractable text.** In the text
layer they render as bare "(1)", "(2)", "(3)" with no formula. Their algebraic content is
**unreadable** from the PDF. What is recoverable is the paper's own verbal statement of the limiting
behaviour (ms p. 20 / PDF p. 21, verbatim):
> "When the term B[RSH)]T < 1 (see equation 1) the reaction approximates to pseudo first-order
> behavior with a linear dependence between kobs and the thiol concentrations. … In such situations,
> A will be equal to an apparent second order rate constant."

Scheme 1, Scheme 2 and Figure 10 are structure drawings with no extractable numbers.

## 3. TABLE 1 — proteins

**Anchor: Table 1, manuscript p. 14 (PDF page 15).** Title as printed:
"Table 1 Apparent second-order rate constant (k2) for reaction of 4MBQ with HSA, BSA, NEM-BSA, and
α-lactalbumin at 25 °C, thiol concentration determined for each protein, and the second-order rate
constants (k2,SH) corrected for reaction with Cys residues (n.d.; not determined)."

Column headers as printed: `Protein` | `Thiol concentration (mol/mol protein)` | `pH` | `k2 (M-1 s-1)`
Units quoted exactly as printed: `M-1 s-1` (i.e. M⁻¹ s⁻¹); thiol concentration in `mol/mol protein`.

| Protein | Thiol concentration (mol/mol protein) | pH | k2 (M⁻¹ s⁻¹) | prov. |
|---|---|---|---|---|
| HSA | 0.21 ± 0.00 | 7.0 | (4.8 ± 0.2) x 10³ | [M] |
| HSA (same row, second line) | — | 7.0 | k2,SH (2.3 ± 0.1) x 10⁴ | [F] (author-computed correction, see §3.1) |
| BSA | 0.38 ± 0.00 | 7.0 | (3.1 ± 0.2) x 10⁴ | [M] |
| BSA (same row, second line) | — | 7.0 | k2,SH (8.3 ± 0.5) x 10⁴ | [F] |
| NEM-BSA | 0.16 ± 0.01 | 7.0 | (1.0 ± 0.1) x 10⁴ | [M] |
| α-lactalbumin | n.d. | 7.0 | (4.0 ± 0.2) x 10² | [M] |
| BSA | 0.38 ± 0.00 | 6.5 | (7.7 ± 0.2) x 10³ | [M] |

Raster-verified against PDF p. 15 cell by cell; no cell unreadable.

### 3.1 Author-stated derivation of k2,SH, and my check

Text (ms p. 13 / PDF p. 14), verbatim: "Such experiments gave values of 0.21 mol/mol protein for HSA
and 0.38 mol/mol protein for BSA (Table 1). Based on these concentrations, corrected second-order
rate constants for reaction with the Cys residues, k2,SH, at pH 7.0 were calculated for BSA and HSA
as (8.3 ± 0.5) × 10⁴ and (2.3 ± 0.1) × 10⁴ M⁻¹s⁻¹, respectively."

| check | arithmetic | printed | verdict |
|---|---|---|---|
| BSA k2,SH | 3.1 × 10⁴ / 0.38 = 8.16 × 10⁴ | 8.3 × 10⁴ | reproduces to 2% [Z] |
| HSA k2,SH | 4.8 × 10³ / 0.21 = 2.29 × 10⁴ | 2.3 × 10⁴ | reproduces exactly [Z] |
| NEM-BSA k2,SH (**never printed by the authors**) | 1.0 × 10⁴ / 0.16 = **6.3 × 10⁴ M⁻¹ s⁻¹** | — | [Z] — do not cite as measured |

### 3.2 INTERNAL INCONSISTENCY — the NEM-blocking percentage

The **abstract** (ms p. 1 / PDF p. 2) says: "Reaction of Cys-34 of BSA with N-ethylmaleimide reduced
the thiol concentration by ~59%, which resulted in a decrease in k2 by a similar percentage,
consistent with rapid adduction at Cys-34."

The **body** (ms p. 14 / PDF p. 15, raster-verified) says: "The thiol concentration of the NEM-BSA
preparation was found to be 0.16 ± 0.01 mol/mol protein, corresponding to a thiol-blocking of ~59%.
The k2 of the reaction of 4MBQ with NEM-BSA at pH 7.0 was decreased by a similar percentage (~67%)
compared to non-blocked."

My arithmetic [Z]: thiol drop = (0.38 − 0.16)/0.38 = **57.9%** (printed "~59%"); k2 drop =
(3.1 × 10⁴ − 1.0 × 10⁴)/3.1 × 10⁴ = **67.7%** (printed "~67%"). The two percentages are not equal
(58% vs 68%); the abstract's "a similar percentage" without a number, and the body's "~67%", are
both quoted rather than reconciled. Both are reported here as printed.

## 4. TABLE 2 — low-molecular-mass thiols, tetrahedral-intermediate parameters

**Anchor: Table 2, manuscript p. 20 (PDF page 21).** Title as printed:
"Table 2 Equilibrium constant (K) and rate constant (k) for reaction of 4MBQ with Nα-acetyl-L-Cys,
L-Cys and GSH. K401 nm and k401 nm were measured by following the loss of 4MBQ at 401 nm, K294 nm and
k294 nm were measured by following the adduct formation at 294 nm."

Column headers as printed: `Thiols` | `K401 nm` | `k401 nm (s-1)` | `K294 nm` | `k294 nm(s-1)`
Units as printed: K columns carry **no units at all**; k columns carry `s-1`.

| Thiols | K401 nm | k401 nm (s⁻¹) | K294 nm | k294 nm (s⁻¹) | prov. |
|---|---|---|---|---|---|
| Nα-Acetyl-L-Cys | (1.6 ± 0.1)×10⁻⁴ | (1.0 ± 0.7)×10³ | (1.8 ± 0.2)×10⁻⁴ | (1.2 ± 0.2)×10³ | [F] |
| L-Cys | (3.3 ± 0.9)×10⁻⁴ | (2.2 ± 0.8)×10³ | (3.2 ± 0.4)×10⁻⁴ | (2.7 ± 0.4)×10³ | [F] |
| GSH | (4.2 ± 1.1)×10⁻⁴ | (1.3 ± 0.4)×10³ | (5.9 ± 0.4)×10⁻⁴ | (1.3 ± 1.0)×10³ | [F] |

Raster-verified against PDF p. 21. Provenance is **[F]** throughout: these are not measured
quantities but parameters obtained from the A and B fit coefficients of equation 1 via the
(image-only, hence unreadable) equations 2 and 3, using the corresponding pH values.

**Units of K are not printed.** The paper's text calls K "the equilibrium constant … for the
formation of the tetrahedral intermediates" and equations 2/3 are images. K therefore cannot be
assigned a unit from this file; recorded as **unit unstated**. Do NOT assume M⁻¹.

Text corroboration (ms p. 19 / PDF p. 20), verbatim: "The equilibrium constants K for the formation
of the tetrahedral intermediates with 4MBQ varied from 1.6 x 10-4 for Nα-acetyl-L-Cys to 4.2 x 10-4
for GSH (Table 2). The corresponding values of k for decay of the intermediate into products varied
from 1.0 x 103 s-1 for Nα-acetyl-L-Cys to 2.2 x 103 s-1 for L-Cys."

## 5. LOW-MOLECULAR-MASS THIOL SECOND-ORDER RATE CONSTANTS (the "A" values) — main-text only, SI table absent

**Anchor: main text, manuscript p. 20 (PDF page 21), quoting "Tables 2, 4, and 6, supplementary
material" — those SI tables are NOT in this PDF.** Verbatim, raster-verified:
> "The high values of A for Nα-acetyl-L-Cys, L-Cys and GSH (5.2 x 10⁵ M⁻¹ s⁻¹ (pH 6.5), 7.0 x 10⁵
> M⁻¹ s⁻¹ (pH 6.0) and 5.4 x 10⁵ M⁻¹ s⁻¹ (pH 6.0), respectively. (Tables 2, 4, and 6, supplementary
> material) indicate that these low-molecular-mass thiols react overall much faster (~68 fold or
> higher) with 4MBQ than the thiol containing proteins, where k2 for BSA at pH 6.5 was 7.7 x 10³
> M⁻¹ s⁻¹ (Table 1)."

| nucleophile | A (= apparent 2nd-order k2 in the linear limit) | units as printed | pH | T | prov. |
|---|---|---|---|---|---|
| Nα-acetyl-L-Cys | 5.2 × 10⁵ | M⁻¹ s⁻¹ | 6.5 | 25 °C | [F], SI Table 2 (absent) |
| L-Cys | 7.0 × 10⁵ | M⁻¹ s⁻¹ | 6.0 | 25 °C | [F], SI Table 4 (absent) |
| GSH | 5.4 × 10⁵ | M⁻¹ s⁻¹ | 6.0 | 25 °C | [F], SI Table 6 (absent) |

No error bars are quoted for these three A values in the main text. **The pH values differ between
the three thiols (6.5 vs 6.0 vs 6.0) — they are NOT an isothermal, iso-pH comparison series.**

My check [Z]: 5.2 × 10⁵ / 7.7 × 10³ = **67.5**, reproducing the paper's "~68 fold".

Discussion restatement (ms p. 29 / PDF p. 30), verbatim: "At pH 7.0, k2 for Nα-acetyl-L-Cys was too
fast to be determined accurately, but at pH 6.5 k2 (corresponding to the A value; Table 2,
supplementary material) was 5.2 × 10⁵ M⁻¹s⁻¹ and since this value is expected to be higher at pH
7.0, the apparent second-order rate constant is assumed to be at least 5 × 10⁵-fold greater than
that for reaction with Nα-acetyl-L-Lys (k2 0.9 ± 0.1 M⁻¹s⁻¹)."

**Note the cross-pH comparison the authors make here:** the "5 × 10⁵-fold" figure divides a pH-6.5
thiol A value by a pH-7.0 amine k2. My check [Z]: 5.2 × 10⁵ / 0.9 = 5.8 × 10⁵. Recorded as the
authors' own construction, not as a measured ratio at a single pH.

## 6. TABLE 3 — amines (Lys and its N-acetyl derivatives)

**Anchor: Table 3, manuscript p. 22 (PDF page 23).** Title as printed:
"Table 3 Second-order rate constants (k2) for the reaction of 4MBQ with L-Lys, Nε-acetyl-L-Lys, and
Nα-acetyl-L-Lys at 25 °C. k2^401nm was determined by following the loss of 4MBQ at 401 nm; k2^adduct
was determined by following the formation of adduct at 499 nm for Nα-acetyl-L-Lys, and at 508 nm for
Nε-acetyl-L-Lys and L-Lys."

Column headers as printed (two header rows):
row 1 `Amines` | `pH 7.0` | `pH 8.0` | `pH 8.0`
row 2 (blank) | `k2^401nm (M-1 s-1)` | `k2^401nm (M-1 s-1)` | `k2^adduct (M-1 s-1)`

| Amines | pH 7.0, k2^401nm (M⁻¹ s⁻¹) | pH 8.0, k2^401nm (M⁻¹ s⁻¹) | pH 8.0, k2^adduct (M⁻¹ s⁻¹) | prov. |
|---|---|---|---|---|
| L-Lys | 8.4 ± 0.5 | 45.7 ± 3.5 | 16.3 ± 0.9 | [M] |
| Nε-acetyl-L-Lys | 4.0 ± 0.4 | 25.8 ± 0.8 | 9.8 ± 0.8 | [M] |
| Nα-acetyl-L-Lys | 0.9 ± 0.1 | 24.8 ± 2.0 | 4.7 ± 0.5 | [M] |

Raster-verified against PDF p. 23; no cell unreadable. Reactivity order stated in text: "Lys >
Nε-acetyl-L-Lys > Nα-acetyl-L-Lys" at pH 7.0 — which the table supports (8.4 > 4.0 > 0.9). At pH 8.0
the two acetyl derivatives are statistically indistinguishable in the 401 nm column (25.8 ± 0.8 vs
24.8 ± 2.0) and the stated ordering no longer separates them; the paper does not comment on this.

[Z] pH-sensitivity ratios, my computation, never printed by the authors:
L-Lys 45.7/8.4 = **5.4×** per pH unit; Nε-acetyl-L-Lys 25.8/4.0 = **6.5×**; Nα-acetyl-L-Lys
24.8/0.9 = **27.6×** over pH 7.0 → 8.0.

## 7. OTHER NITROGEN-NUCLEOPHILE RATE CONSTANTS — main text only (SI Tables 8, 9 absent)

**Anchor: main text, manuscript pp. 20–21 (PDF pp. 21–22).** Verbatim:
> "The reaction of 4MBQ with the guanidine group of Nα-acetyl-L-Arg was too slow to be determined
> accurately under the conditions employed. At pH 6.5, kobs for the decay of 4MBQ (0.04 mM) in the
> presence of 50 mM Nα-acetyl-L-Arg was (7.5 ± 0.2) × 10-2 s-1, which was not significantly different
> from the value of kobs for spontaneous decay of 4MBQ ((7.9 ± 0.2) × 10-2 s-1, p > 0.05). Similar data
> were obtained at pH 7.0 (data not shown)."
> "The apparent second-order rate constant, k2, for loss of 4MBQ in the presence of Gly at pH 6.5 was
> 0.7 ± 0.1 M-1 s-1 (Table 8, supplementary material), which is 7 × 105 times lower than for Nα-acetyl-
> L-Cys (5.2 x 105 M-1 s-1) at the same pH. At pH 7.0, k2 for Gly increased to 2.0 ± 0.3 M-1 s-1. The
> reaction of 4MBQ with L-Lys at pH 7.0 (k2 8.4 ± 0.5 M-1 s-1) was faster than for Gly, but still
> significantly slower than that for any of the thiols (Table 9, supplementary material)."

| nucleophile | quantity | value | units as printed | pH | T | prov. |
|---|---|---|---|---|---|---|
| Nα-acetyl-L-Arg (50 mM) | kobs (4MBQ decay) | (7.5 ± 0.2) × 10⁻² | s⁻¹ | 6.5 | 25 °C | [M] — **upper bound only**, not distinguishable from blank |
| (none — blank) | kobs (spontaneous 4MBQ decay) | (7.9 ± 0.2) × 10⁻² | s⁻¹ | 6.5 | 25 °C | [M] |
| Gly | k2 | 0.7 ± 0.1 | M⁻¹ s⁻¹ | 6.5 | 25 °C | [F], SI Table 8 (absent) |
| Gly | k2 | 2.0 ± 0.3 | M⁻¹ s⁻¹ | 7.0 | 25 °C | [F], SI Table 9 region (absent) |
| L-Lys | k2 | 8.4 ± 0.5 | M⁻¹ s⁻¹ | 7.0 | 25 °C | [M] (= Table 3) |

**No k2 can be assigned to Nα-acetyl-L-Arg.** The measured kobs with 50 mM Arg is *lower* than the
blank, so the reaction is below the detection floor of the method. Any Arg rate constant must be
recorded as "< blank-limited", never as a number.

My check [Z]: 5.2 × 10⁵ / 0.7 = 7.4 × 10⁵, reproducing the paper's "7 × 10⁵ times lower".
Also: the non-zero intercept in Figure 3 is explicitly attributed to this spontaneous decay ("The
non-zero intercept on the vertical axis is due to the slow spontaneous decay of 4MBQ determined in
the absence of added protein", Fig. 3 caption, ms p. 11 / PDF p. 12).

## 8. THIOL pKa VALUES USED BY THE AUTHORS — all cited, none measured

**Anchor: Discussion, manuscript p. 29 (PDF page 30).** Verbatim: "…which may be explained, at least
in part, by the high pKa of Nα-acetyl-L-Cys (thiol pKa 9.5) [45]. However, the values of K determined
for L-Cys (thiol pKa 8.3) [45] were smaller than those of GSH (thiol pKa 8.8) [45]."

| compound | thiol pKa | prov. |
|---|---|---|
| Nα-acetyl-L-Cys | 9.5 | **[C]** — Winterbourn & Metodiewa, *Free Radic. Biol. Med.* **27** (1999) 322-328 (ref. 45) |
| L-Cys | 8.3 | **[C]** — same source |
| GSH | 8.8 | **[C]** — same source |

The authors then note the correlation **fails**: "the values of K determined for L-Cys (thiol pKa 8.3)
were smaller than those of GSH (thiol pKa 8.8)". No structure–reactivity regression is fitted
anywhere in this paper.

## 9. TABLE 4 — LC-ESI-MS adduct identification (no kinetics)

**Anchor: Table 4, manuscript p. 25 (PDF page 26).** Title as printed:
"Table 4 ESI (positive ion mode) MS data for thiol-phenol and amine-quinone adducts."
Column headers as printed: `Compound` | `Molecular Mass` | `Retention time (min)` | `m/z (relative intensity, %)` | `Proposed ion structure`

| Compound | Molecular Mass | Retention time (min) | m/z (relative intensity, %) | Proposed ion structure | prov. |
|---|---|---|---|---|---|
| L-Cys-phenol adduct | 243.2 | 12.7 | 244.0 | [M+H]⁺ | [M] |
| Nα-acetyl-L-Cys-phenol adduct | 285.3 | 5.5 | 286.0 | [M+H]⁺ | [M] |
| GSH-phenol adduct | 429.4 | 6.8 | 430.1 | [M+H]⁺ | [M] |
| Gly-quinone adduct | 195.2 | 6.3 | 196.0 (100) | [M+H]⁺ | [M] |
| Gly-quinone adduct (2nd line) | — | — | 218.0 (10) | [M+Na]⁺ | [M] |
| L-Lys-quinone adduct | 266.3 | 7.2 | 267.1 | [M+H]⁺ | [M] |
| Nα-acetyl-L-Lys-quinone adduct | 308.3 | 9.7 | 309.1 (100) | [M+H]⁺ | [M] |
| Nα-acetyl-L-Lys-quinone adduct (2nd line) | — | — | 331.1 (15) | [M+Na]⁺ | [M] |
| Nε-acetyl-L-Lys-quinone adduct | 308.3 | 7.9 | 309.1 (100) | [M+H]⁺ | [M] |
| Nε-acetyl-L-Lys-quinone adduct (2nd line) | — | — | 331.1 (58) | [M+Na]⁺ | [M] |

## 10. TABLE 5 — MS/MS fragmentation (no kinetics)

**Anchor: Table 5, manuscript p. 26 (PDF page 27).** Title as printed:
"Table 5 ESI MS/MS (positive ion mode) data for thiol-phenol and amine-quinone adduct."
Column headers as printed: `Compound` | `Precursor ion m/z` | `Product ions, m/z (relative intensity, %)`

| Compound | Precursor ion m/z | Product ions, m/z (relative intensity, %) | prov. |
|---|---|---|---|
| L-Cys-phenol adduct | 244.0 | 226.9 (100), 155.0 (53) | [M] |
| Nα-acetyl-L-Cys- phenol adduct | 286.0 | 268.0 (100), 244.0 (87), 227.0 (19), 198.0 (16), 181.0 (64), 161.9 (99) | [M] |
| GSH-phenol adduct | 430.1 | 355.0 (4), 301.0 (100), 284.0 (46), 130.0 (34) | [M] |
| Gly-quinone adduct | 196.0 | 150.0 (100), 122.0 (20) | [M] |
| L-Lys-quinone adduct | 267.1 | 250.0 (6), 205.8 (4), 128.0 (100) | [M] |
| Nα-acetyl-L-Lys-quinone adduct | 309.1 | 204.1 (6), 170.0 (10), 128.0 (100) | [M] |
| Nε-acetyl-L-Lys-quinone adduct | 309.1 | 265.0 (8), 170.0 (4), 128.0 (100) | [M] |

## 11. FIGURE-ONLY QUANTITIES

| figure | anchor | content | usable numbers |
|---|---|---|---|
| Figure 1 | ms p. 9 / PDF p. 10 | UV-Vis of 1 mM 4MC, 4MBQ, and re-reduced 4MBQ, 22 °C, pH 4.5 | λmax 280 nm (4MC), 401 nm (4MBQ); "around 85 % conversion back to 4MC" [M, text] |
| Figure 2 | ms p. 10 / PDF p. 11 | 0.2 mM HSA + 0.04 mM 4MBQ, pH 7.0, 0.2 M phosphate, 25.0 °C, 10 s | no tabulated values; [fig] only |
| Figure 3 | ms p. 11 / PDF p. 12 | kobs vs [HSA]T, 0.2–0.5 mM | slope = k2 4.8 × 10³ M⁻¹s⁻¹ (already in Table 1); non-zero intercept ≈ spontaneous decay |
| Figure 4 | ms p. 12 / PDF p. 13 | 4MBQ decay vs BSA / NEM-BSA / α-lactalbumin, all 0.28 mM, pH 7.0 | no tabulated values |
| Figure 5 | ms p. 15–16 / PDF p. 16–17 | Nα-acetyl-L-Cys 0.5 mM + 4MBQ 0.04 mM, pH **5.0**, 25.0 °C | adduct λmax 294 nm |
| Figure 6 | ms p. 18 / PDF p. 19 | kobs vs [Nα-acetyl-L-Cys]T at **pH 5.0**, x-axis to 6 × 10⁻³ M, y-axis to 80 s⁻¹ (401 nm) and 100 s⁻¹ (294 nm) | saturating (non-linear) shape; A and B values themselves are in the absent SI |
| Figure 7 | ms p. 19 / PDF p. 20 | A and B vs 1/[H⁺], x-axis to 3 × 10⁶ M⁻¹; A-axis to 6–8 × 10⁵, B-axis to 2 × 10³ | axis ranges only; the fitted A(pH), B(pH) numbers are in the absent SI |
| Figure 8 | ms p. 23 / PDF p. 24 | L-Lys 2.5 mM + 4MBQ 0.04 mM, pH **8.0**, 20 s; insert shows a lag of ~0.5 s | lag phase < 0.5 s [fig] |
| Figure 9 | ms p. 24 / PDF p. 25 | amine-quinone adduct spectra, 25 mM amine + 0.2 mM 4MBQ, pH 7.0, 22 °C | λmax ~499 nm (Nα-acetyl-L-Lys); 504–508 nm (Gly, Nε-acetyl-L-Lys) [M, text] |
| Figure 11 | ms p. 30 / PDF p. 31 | 4MBQ 0.04 mM vs BSA 0.5 mM, Nα-acetyl-L-Cys 0.5 mM, Nα-acetyl-L-Arg 0.5 mM, Gly 0.5 mM, pH **6.5**, 25 °C | qualitative ordering only |

I have **not** digitised any figure point in this dossier; nothing in the parameter table below is
tagged [fig].

## 12. DEFECTS AND CAVEATS TO CARRY FORWARD

1. **No temperature series → no Ea from this paper.** Everything is 25 °C. Any Arrhenius transfer to
   food temperatures is unlicensed by this source.
2. **The supplementary tables (SI Tables 2–9) are not in this PDF.** The A/B fit coefficients per pH,
   the kobs raw data, and the Gly/Lys concentration series exist only as the numbers the main text
   quotes. Any wave that wants the full A(pH) curve must retrieve the SI separately.
3. **Equations 1, 2, 3 are images with no text layer — algebraically unreadable from this file.**
   The K and k values in Table 2 therefore cannot be re-derived or unit-checked here.
4. **K in Table 2 has no printed units.** Do not assume M⁻¹.
5. **± values are undefined** (SD? SE? fit error?) and **n is never stated**.
6. **The stated pH range "4.5 to 7.0" (ms p. 9) contradicts Table 3, which reports pH 8.0.**
7. The three LMM-thiol A values are at **three different pH values** (6.5, 6.0, 6.0) — not a
   like-for-like thiol series.
8. Nα-acetyl-L-Arg gives **no rate constant**; its kobs is below the spontaneous-decay blank.
9. The thiol pKa values are **[C] cited from ref. 45**, not measured here, and the paper itself
   reports that K does not track pKa monotonically (L-Cys pKa 8.3 gives smaller K than GSH pKa 8.8).
10. Abstract "~59%" vs body "~67%" for the NEM-induced k2 drop (see §3.2).

## NEW-PARAMETER TABLE (consolidated)

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| k2, 4MBQ + HSA | (4.8 ± 0.2) × 10³ | M⁻¹ s⁻¹ | 25 °C, 0.2 M phosphate, pH 7.0, [4MBQ] 0.04 mM, [HSA] 0.2–0.5 mM | Table 1, ms p. 14 / PDF p. 15 | [M] |
| k2,SH, 4MBQ + HSA Cys-34 | (2.3 ± 0.1) × 10⁴ | M⁻¹ s⁻¹ | as above, corrected by 0.21 mol SH/mol protein | Table 1, ms p. 14 / PDF p. 15 | [F] |
| k2, 4MBQ + BSA | (3.1 ± 0.2) × 10⁴ | M⁻¹ s⁻¹ | 25 °C, 0.2 M phosphate, pH 7.0 | Table 1, ms p. 14 / PDF p. 15 | [M] |
| k2,SH, 4MBQ + BSA Cys-34 | (8.3 ± 0.5) × 10⁴ | M⁻¹ s⁻¹ | as above, corrected by 0.38 mol SH/mol protein | Table 1, ms p. 14 / PDF p. 15 | [F] |
| k2, 4MBQ + NEM-BSA | (1.0 ± 0.1) × 10⁴ | M⁻¹ s⁻¹ | 25 °C, pH 7.0, 0.16 mol SH/mol protein (~59% blocked) | Table 1, ms p. 14 / PDF p. 15 | [M] |
| k2,SH, 4MBQ + NEM-BSA | 6.3 × 10⁴ | M⁻¹ s⁻¹ | 1.0 × 10⁴ / 0.16 | Table 1, ms p. 14 / PDF p. 15 | **[Z]** — never printed |
| k2, 4MBQ + α-lactalbumin | (4.0 ± 0.2) × 10² | M⁻¹ s⁻¹ | 25 °C, pH 7.0, no free Cys | Table 1, ms p. 14 / PDF p. 15 | [M] |
| k2, 4MBQ + BSA | (7.7 ± 0.2) × 10³ | M⁻¹ s⁻¹ | 25 °C, 0.2 M phosphate, **pH 6.5** | Table 1, ms p. 14 / PDF p. 15 | [M] |
| BSA thiol content | 0.38 ± 0.00 | mol/mol protein | DTNB, Tris pH 8.0, 25 °C, 30 min | Table 1, ms p. 14 / PDF p. 15 | [M] |
| HSA thiol content | 0.21 ± 0.00 | mol/mol protein | as above | Table 1, ms p. 14 / PDF p. 15 | [M] |
| NEM-BSA thiol content | 0.16 ± 0.01 | mol/mol protein | as above | Table 1, ms p. 14 / PDF p. 15 | [M] |
| K401 nm, 4MBQ + Nα-acetyl-L-Cys | (1.6 ± 0.1) × 10⁻⁴ | *(no units printed)* | 25 °C, 0.2 M phosphate, pH-resolved fit | Table 2, ms p. 20 / PDF p. 21 | [F] |
| k401 nm, Nα-acetyl-L-Cys intermediate decay | (1.0 ± 0.7) × 10³ | s⁻¹ | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| K294 nm, 4MBQ + Nα-acetyl-L-Cys | (1.8 ± 0.2) × 10⁻⁴ | *(no units printed)* | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| k294 nm, Nα-acetyl-L-Cys | (1.2 ± 0.2) × 10³ | s⁻¹ | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| K401 nm, 4MBQ + L-Cys | (3.3 ± 0.9) × 10⁻⁴ | *(no units printed)* | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| k401 nm, L-Cys | (2.2 ± 0.8) × 10³ | s⁻¹ | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| K294 nm, 4MBQ + L-Cys | (3.2 ± 0.4) × 10⁻⁴ | *(no units printed)* | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| k294 nm, L-Cys | (2.7 ± 0.4) × 10³ | s⁻¹ | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| K401 nm, 4MBQ + GSH | (4.2 ± 1.1) × 10⁻⁴ | *(no units printed)* | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| k401 nm, GSH | (1.3 ± 0.4) × 10³ | s⁻¹ | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| K294 nm, 4MBQ + GSH | (5.9 ± 0.4) × 10⁻⁴ | *(no units printed)* | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| k294 nm, GSH | (1.3 ± 1.0) × 10³ | s⁻¹ | as above | Table 2, ms p. 20 / PDF p. 21 | [F] |
| A (apparent k2), 4MBQ + Nα-acetyl-L-Cys | 5.2 × 10⁵ | M⁻¹ s⁻¹ | 25 °C, **pH 6.5** | main text ms p. 20 / PDF p. 21 (SI Table 2, absent) | [F] |
| A (apparent k2), 4MBQ + L-Cys | 7.0 × 10⁵ | M⁻¹ s⁻¹ | 25 °C, **pH 6.0** | main text ms p. 20 / PDF p. 21 (SI Table 4, absent) | [F] |
| A (apparent k2), 4MBQ + GSH | 5.4 × 10⁵ | M⁻¹ s⁻¹ | 25 °C, **pH 6.0** | main text ms p. 20 / PDF p. 21 (SI Table 6, absent) | [F] |
| k2, 4MBQ + L-Lys | 8.4 ± 0.5 | M⁻¹ s⁻¹ | 25 °C, pH 7.0, 401 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2, 4MBQ + L-Lys | 45.7 ± 3.5 | M⁻¹ s⁻¹ | 25 °C, pH 8.0, 401 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2^adduct, 4MBQ + L-Lys | 16.3 ± 0.9 | M⁻¹ s⁻¹ | 25 °C, pH 8.0, 508 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2, 4MBQ + Nε-acetyl-L-Lys | 4.0 ± 0.4 | M⁻¹ s⁻¹ | 25 °C, pH 7.0, 401 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2, 4MBQ + Nε-acetyl-L-Lys | 25.8 ± 0.8 | M⁻¹ s⁻¹ | 25 °C, pH 8.0, 401 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2^adduct, 4MBQ + Nε-acetyl-L-Lys | 9.8 ± 0.8 | M⁻¹ s⁻¹ | 25 °C, pH 8.0, 508 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2, 4MBQ + Nα-acetyl-L-Lys | 0.9 ± 0.1 | M⁻¹ s⁻¹ | 25 °C, pH 7.0, 401 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2, 4MBQ + Nα-acetyl-L-Lys | 24.8 ± 2.0 | M⁻¹ s⁻¹ | 25 °C, pH 8.0, 401 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2^adduct, 4MBQ + Nα-acetyl-L-Lys | 4.7 ± 0.5 | M⁻¹ s⁻¹ | 25 °C, pH 8.0, 499 nm | Table 3, ms p. 22 / PDF p. 23 | [M] |
| k2, 4MBQ + Gly | 0.7 ± 0.1 | M⁻¹ s⁻¹ | 25 °C, pH 6.5 | main text ms p. 21 / PDF p. 22 (SI Table 8, absent) | [F] |
| k2, 4MBQ + Gly | 2.0 ± 0.3 | M⁻¹ s⁻¹ | 25 °C, pH 7.0 | main text ms p. 21 / PDF p. 22 | [F] |
| kobs, 4MBQ + Nα-acetyl-L-Arg | (7.5 ± 0.2) × 10⁻² | s⁻¹ | 25 °C, pH 6.5, 50 mM Arg, 0.04 mM 4MBQ — **≤ blank, no k2 derivable** | main text ms p. 20 / PDF p. 21 | [M] |
| kobs, spontaneous 4MBQ decay | (7.9 ± 0.2) × 10⁻² | s⁻¹ | 25 °C, pH 6.5, no nucleophile | main text ms p. 20 / PDF p. 21 | [M] |
| ε395, 4MBQ | 1350 | M⁻¹ cm⁻¹ | aqueous | Methods ms p. 5 / PDF p. 6 | [C] (ref. 28, Whitaker et al.) |
| thiol pKa, Nα-acetyl-L-Cys | 9.5 | (dimensionless) | — | Discussion ms p. 29 / PDF p. 30 | [C] (ref. 45) |
| thiol pKa, L-Cys | 8.3 | (dimensionless) | — | Discussion ms p. 29 / PDF p. 30 | [C] (ref. 45) |
| thiol pKa, GSH | 8.8 | (dimensionless) | — | Discussion ms p. 29 / PDF p. 30 | [C] (ref. 45) |

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| Table 1, **pH 7.0 rows only** (HSA 4.8 × 10³; BSA 3.1 × 10⁴; NEM-BSA 1.0 × 10⁴; α-lactalbumin 4.0 × 10²) + the DTNB thiol contents 0.21 / 0.38 / 0.16 | **FIT** | pH (arm A) | This is the only self-consistent iso-pH, iso-T quartet in the paper. It fixes the thiol-vs-non-thiol reactivity ratio (BSA/α-lactalbumin = 78×) and the per-thiol normalisation, which is exactly the quinone–thiol branch parameter the trunk lacks. The NEM row is an internal control on the same axis, not an independent datum. |
| Table 1, **BSA pH 6.5 row (7.7 × 10³ M⁻¹ s⁻¹)** | **HOLD-OUT** | pH (arm B) | Same protein, same T, same buffer, only pH moves. Once the pH-7.0 arm is fit, predicting 7.7 × 10³ at pH 6.5 is a genuine one-variable extrapolation of the thiolate-fraction term. Cheap, clean, and orthogonal. |
| Table 3, **pH 7.0 column** (L-Lys 8.4; Nε-acetyl-L-Lys 4.0; Nα-acetyl-L-Lys 0.9) | **FIT** | pH (arm A) | Sets the amine branch at the same pH as the FIT protein arm, so the thiol:amine selectivity is fixed by data taken at one pH rather than by the paper's own cross-pH "5 × 10⁵-fold" construction. |
| Table 3, **pH 8.0 columns** (both k2^401nm and k2^adduct; 6 numbers) | **HOLD-OUT** | pH (arm B) | Genuine extrapolation on the amine branch, and it additionally tests whether the model reproduces the *ordering collapse* between the two acetyl derivatives at pH 8.0 (25.8 vs 24.8), which a naive pKa-shift model will get wrong. |
| Gly k2: 0.7 (pH 6.5) and 2.0 (pH 7.0) | **HOLD-OUT** (both) | nucleophile identity + pH | Gly was never used to set anything; it is the smallest amine and the largest lever arm on the amine branch. Holding both out gives a two-point pH slope check on a nucleophile absent from the FIT set. |
| Table 2 (K and k for the tetrahedral intermediate, all 12 numbers) | **DO NOT USE** | — | Units of K are not printed and equations 1–3 are image-only and unreadable in this PDF. These are fit coefficients of a mechanism the file does not let us reconstruct. They must not enter any parameter until the SI is retrieved. |
| A values (5.2 / 7.0 / 5.4 × 10⁵ M⁻¹ s⁻¹) | **HOLD-OUT, low weight** | pH is confounded | Each is at a different pH (6.5 / 6.0 / 6.0), so they are not a comparable series. Usable only as an order-of-magnitude ceiling on LMM-thiol reactivity, never as a three-point structure–reactivity fit. |
| Nα-acetyl-L-Arg | **NOT A DATUM** | — | kobs is below the spontaneous-decay blank. Record as "Arg branch negligible", never as a rate constant. |
| Thiol pKa 9.5 / 8.3 / 8.8 | **not fittable here** | — | **[C]**, from Winterbourn & Metodiewa 1999. If the trunk already carries thiol pKa values from another source, these must not be double-counted as independent evidence. |

**Circularity risks flagged.**
(i) The k2,SH rows in Table 1 are **not independent measurements** — they are k2 divided by the DTNB
thiol content in the same row. Fitting both k2 and k2,SH from Table 1 would double-count one
measurement. Use k2 + thiol content, or k2,SH alone, never all three.
(ii) The NEM-BSA row is derived from the same BSA preparation as the BSA row; treating it as an
independent protein would inflate the effective n.
(iii) The "5 × 10⁵-fold thiol:amine" headline mixes a pH-6.5 thiol value with a pH-7.0 amine value.
If the model is scored against that ratio it will be scored against an artefact of the authors'
comparison, not against a measurement.
(iv) **No Ea is available from this paper at all.** Any Arrhenius extrapolation of these constants to
food temperatures is unsupported by this source and must be declared as such.

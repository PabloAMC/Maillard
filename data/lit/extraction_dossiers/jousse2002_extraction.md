# Jousse, Jongen, Agterof, Russell & Braat 2002 (10.1111/j.1365-2621.2002.tb08772.x) — per-paper extraction 2026-09-07

**Source PDF:** `data/articles/jousse2002.pdf` (9 pp., JFS 67(7):2534–2542). Born-digital (PageMaker 6.5 →
Acrobat Distiller 4.05, created 24 Sep 2002); text layer complete and reliable for prose, **but it flattens
every superscript** (`7.4 × 1014` = 7.4 × 10¹⁴, `10-7` = 10⁻⁷) and drops Eq. 1 and Eq. 2 entirely.
Read method: **both** — full text layer, **plus** 200 dpi rasters of pp. 2535–2541 (`j_p2..j_p8`) for
Figure 1 (the scheme), Eq. 1 (the ODE system), Eq. 2 (heating lag), Tables 1, 2, 4, 5 and the Figure 8
caption. **Every numeric cell below was read from the raster, not the text layer.**

### ★ HEADLINE: this is the corpus's only *lumped-class* mass-action network for Maillard **volatiles** (PY / FU / C / PZ) with a full (prefactor, Ea) set — but the set is a **hand-tuned Excel fit in a glycerol matrix**, and the literature-derived prefactors are printed **inconsistently between Tables 1, 2 and 5** (§0.2). Use the Ea values as priors; do not transplant the prefactors.

---

## §0. IDENTITY

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/jousse2002.pdf` — **584,239 bytes**, 9 pages, PDF 1.3, 585 × 783 pt | `ls`, `pdfinfo` |
| SHA-256 | `e80432a54178a9964f2fb6a210d504e3b8f2e9a19fb11885674a05c9dc4201e5` | `shasum -a 256` |
| **title** | ***"Simplified Kinetic Scheme of Flavor Formation by the Maillard Reaction"*** | p.2534 raster ✔; PDF metadata title identical |
| **authors** | **F. Jousse, T. Jongen, W. Agterof, S. Russell, and P. Braat** | p.2534 ✔ |
| **affiliation** | **Unilever Research Vlaardingen, Olivier van Noortlaan 120, 3133AT Vlaardingen, The Netherlands** (all five); inquiries to Jousse, `fabien.jousse@unilever.com` | p.2542 |
| **venue** | *Journal of Food Science* **Vol. 67, Nr. 7, 2002, pp. 2534–2542**, section "Food Chemistry and Toxicology" | running headers |
| **DOI** | **10.1111/j.1365-2621.2002.tb08772.x** (given in the brief; **[NEG] no DOI is printed anywhere in the PDF** — 2002 JFS issues carried none) | whole-document sweep |
| manuscript history | **MS 20010573 Submitted 10/11/01, Accepted 12/18/01, Received 12/21/01** | p.2542 |
| typesetter stamp | `jfsv67n7p2534-2542ms20010573-TLS-TO RRD.p65 … 9/24/2002, 2:16 PM` on every page; p.2542 header still carries the proof query *"(IS THIS ACCEPTABLE?)"* | every page |
| **funding** | **European Union Marie Curie Fellowship, contract nr HPMI-CT-1999-00016** "Modeling of transport phenomena (momentum, heat and mass transfer) during food processing" | p.2542 |
| conflicts | **[NEG] none stated** (industrial authorship, Unilever) | — |
| SI / data deposit | **[NEG] none** — no supplementary material; no raw concentration table for the glycerol experiment (Fig. 10–11 are the only data) | — |
| PDF character | born-digital; `Suspects: no`; figures are raster with no numeric text layer | `pdfinfo` |

### §0.1 TEXT-LAYER vs RASTER

Exponents in Tables 1, 2, 5 and the Fig. 8 caption are flattened in the text layer (`1014`,
`10-13`) and Eq. 1 / Eq. 2 are absent from it entirely; all were transcribed from the raster.
Table 4 names are identical in both. Only glyph slips ("a-dicarbonyl", "El+", "70 Ev"); no
digit is affected.

### §0.2 THE PAPER'S OWN INTERNAL INCONSISTENCIES `[Z]`

| # | inconsistency | evidence | severity |
|---|---|---|---|
| **C1** | ★ **R₂ prefactor: Table 1 prints 7.4 × 10¹⁴; Tables 2 and 5 ("No fit") print 7.4 × 10¹³** — one order of magnitude apart for the same regression. Recomputing Table 2's own R(20/100/180 °C) columns from (7.4 × 10¹³, 128.8) gives 8.3 × 10⁻¹⁰ / 6.9 × 10⁻⁵ / 1.05 × 10⁻¹, matching the printed 7.9 × 10⁻¹⁰ / 6.7 × 10⁻⁵ / 1.0 × 10⁻¹ → **Table 2 was computed with 10¹³.** But the Figure 2 regression line reads ≈ 2 × 10⁻⁷ M⁻¹s⁻¹ at 1000/T = 3.2 and ≈ 4 × 10⁻⁴ at 2.68, which (7.4 × 10¹⁴, 128.8) reproduces (2.2 × 10⁻⁷, 6.9 × 10⁻⁴) and 10¹³ does not → **the literature regression is 10¹⁴; the model used 10¹³.** Never stated. | raster ✔ + `[Z]` arithmetic | ★★★ the ARP-formation rate actually simulated is **10× below** the paper's own literature fit |
| **C2** | **R₄ prefactor: Table 2 prints 5.3 × 10³; Table 1 ("ARP decomposition") and Table 5 ("No fit") print 2.3 × 10³.** Table 2's own columns (8.3 × 10⁻⁷ / 8.8 × 10⁻⁵ / 1.8 × 10⁻³) recompute from **2.3 × 10³** (8.6 × 10⁻⁷ / 9.0 × 10⁻⁵ / 1.8 × 10⁻³), not from 5.3 × 10³ (2.0 × 10⁻⁶ / 2.1 × 10⁻⁴ / 4.2 × 10⁻³). | raster ✔ + `[Z]` | ★★ **`5.3 × 10³` is a typo; 2.3 × 10³ is the value used** |
| C3 | Table 5 title says the fit is to *"PY, FU, SA, and PZ"*; the text, Table 4 and Figure 11 fit **PY, FU, C, PZ**. No Strecker aldehyde (SA) was quantified anywhere. | p.2539–2541 | ★ title error |
| C4 | Abstract: *"The scheme comprises 11 reaction steps"*; in practice **R₇ ≡ 0, R₁₀ "Fast", R₁ = 0 in the fit** — the fitted model has **8 numeric rate constants** (R₂–R₆, R₈, R₉, R₁₁), the literature model **8** (R₁–R₆, R₈, R₉). | p.2538, Table 5 | ★ |
| C5 | R₁ (150.0 kJ/mol, 1.7 × 10¹⁴ s⁻¹) appears in Table 2 with **no source, no Figure and no Table 1 row**; it is then set to 0 by the fit. | p.2537, p.2541 | ★ **REFUSE R₁ as a prior** |
| C6 | Table 1 gives ΔR₀, ΔE as "95% confidence interval on the regression"; Table 2/5 carry **no uncertainty at all** on any fitted rate. | — | ★ |

---

## §1. ONE-PARAGRAPH VERDICT — READ BEFORE USING ANY NUMBER HERE

This is a **modelling** paper with a small experiment bolted on. It delivers (i) a 10-species,
11-arrow lumped scheme (Fig. 1 + Eq. 1) whose topology is the closest thing in the corpus to
this repository's trunk lane; (ii) five literature-regression (R₀, Ea) pairs with 95 % CIs
(Table 1) assembled from 66 sources spanning green peas to milk to aqueous glucose–glycine;
(iii) two full parameter sets (Table 2 "no fit" / Table 5 "fit") and (iv) one glucose–alanine
experiment in **glycerol** at 95–140 °C, fitted **"by hand"** in Excel with an **arbitrary ×2**
multiplier on every measured class. **Goodness of fit is never quantified** — no R², no residual,
no CI on any fitted rate; the Fig. 7 "good agreement" is a log–log scatter over **nine orders of
magnitude** and the authors themselves put the sensitivity of the prediction at *"about 1 order
of magnitude."* The usable content is therefore: the **step topology**, the **Ea values with their
CIs as priors** (ARP formation 128.8 ± 29.4; ARP decomposition 52.9 ± 14.0; HMF formation
109.3 ± 31.3; AA loss 84.1 ± 5.8; sugar loss 75.7 ± 10.9 kJ/mol), the **class definitions**
(Table 4), and the authors' own statement of scope: *"The rates given in Table 5 are valid,
strictly speaking, only for a model glucose-alanine system in glycerol at a pH of 6.5."*
Prefactors are matrix-specific and internally inconsistent (C1, C2): **do not transplant them.**
No organic acids, no DMHF-specific, no HMF-specific kinetics exist in the experimental part.

---

## §2. THE SCHEME — verbatim

### §2.1 The four stages (p.2535) `[M]`

> "(1) Condensation of a sugar (S) with an amino acid (AA) to form an Amadori or Heyns
> rearrangement product (ARP or HRP, respectively). Alternatively, the sugar can also directly
> degrade, for example, by a caramelization reaction at high temperature. We postulate that these
> direct degradations are negligible at low to medium temperatures.
> (2) The intermediate ARP or HRP can cyclize to form nitrogen-containing heterocyclic
> compounds, such as pyrroles or pyridines (PY). Alternatively, it may cleave to give rearranged
> sugars (RS), which contain the intact chain of the starting sugar. These RS incorporate the
> 1-desoxy-2,3-diketones, and the 3-desoxy-1,2-diketones, as well as further rearrangement from
> these via ketoenol tautomerization. This cleavage gives back the original amino acid.
> (3) The RS may cyclize into oxygen-containing heterocyclic compounds, such as furans or
> furfurals (FU). It can also break up into α-dicarbonyl fragments (C), which may recombine to
> give FU.
> (4) Dicarbonyl fragments react with the amine group of the amino acid in the Strecker
> degradation mechanism, giving an intermediate (I), common to Strecker aldehydes (SA) and
> pyrazines (PZ). Alternatively, SA can also represent the nitrogen-containing heterocyclic
> compounds such as pyrroline and pyrrolidine coming from the reaction of dicarbonyls with
> proline and hydroxyproline (Mottram 1994).
> In addition, we considered additional 'further reactions' of the flavor components to give
> polymeric brown products. These are modeled as unimolecular reactions with a unique rate
> constant R₁₁ for all molecules."

### §2.2 The 11 steps as drawn in Figure 1 (p.2535, raster) `[M]` — arrow labels as printed

| step | from → to | label on arrow | order (p.2536) | unit |
|---|---|---|---|---|
| R₁ | Sugar (S) → Degradation product (SD) | — | 1st | s⁻¹ |
| R₂ | S + Amino acid (AA) → Amadori rearrangement products (ARP) | "condense" | 2nd | M⁻¹ s⁻¹ |
| R₃ | ARP → Pyrroles (PY) | "cyclize" | 1st | s⁻¹ |
| R₄ | ARP → Rearranged sugars (RS) **+ AA** (regenerated) | "break-up" | 1st | s⁻¹ |
| R₅ | RS → Furans (FU) | "cyclise" | 1st | s⁻¹ |
| R₆ | RS → **2** Carbonyls (C) | "break-up" | 1st | s⁻¹ |
| R₇ | C + C → FU | "condense" | 2nd | M⁻¹ s⁻¹ (**set to 0**) |
| R₈ | C + AA → Strecker intermediate (I) | "Strecker" | 2nd | M⁻¹ s⁻¹ |
| R₉ | I → Strecker aldehydes (SA) | "Strecker" | 1st | s⁻¹ |
| R₁₀ | I + I → Pyrazines (PZ) | "condense" | 2nd | M⁻¹ s⁻¹ (**"Fast"**) |
| R₁₁ | PY, FU, C, SA, PZ → "Melanoidins" | "further reactions" (dashed) | 1st, one constant for all five | s⁻¹ |

> (p.2536) "The rates are given as 'per second' for unimolecular reaction steps R₁, R₃, R₄, R₅,
> R₆, R₉, and R₁₁, and as 'per second per molarity (M)' for bimolecular reactions R₂, R₇, R₈, and
> R₁₀."

### §2.3 Eq. 1 — the ODE system (p.2536, raster; dot = time derivative) `[M]`

```
Ṡ   = −R1·S − R2·S·AA
ȦA  = −R2·S·AA + R4·ARP − R8·C·AA
ȦRP = R2·S·AA − R3·ARP − R4·ARP
ṖY  = R3·ARP − R11·PY
ṘS  = R4·ARP − R5·RS − R6·RS
ḞU  = R5·RS + R7·C·C − R11·FU
Ċ   = 2·R6·RS − R7·C·C − R8·C·AA − R11·C
İ   = R8·C·AA − R9·I − R10·I·I
ṠA  = R9·I − R11·SA
ṖZ  = R10·I·I − R11·PZ
```

`[D]` Notes on the topology: S is consumed only by R₁ and R₂ (no reversibility anywhere); the
amino acid is **regenerated** by R₄ and consumed again by R₈ (a catalytic loop through ARP);
the **factor 2** on R₆ is a C₆ → 2 × C₃ stoichiometry; ARP/HRP are one pool; RS lumps 1- and
3-deoxyosones and their tautomers; SD and "Melanoidins" are terminal sinks with no equation.
[NEG] No water, no pH, no a_w term, no acid species, no reversibility, no sugar isomerisation.

### §2.4 The four volatile classes — Table 4 verbatim (p.2539) `[M]`

Title: *"Composition of the 4 'classes' of compounds used to fit the kinetic scheme. Only
molecules that have been unambiguously identified have been put into each class."* Rt = retention
time (min) on the Figure 9 chromatogram.

| class | Rt (min) | compound as printed |
|---|---|---|
| **Pyrazines (PZ)** | 3.233 | methylpyrazine |
| | 3.435 | 2,5-dimethylpyrazine |
| | 3.522 | 2,6-dimethylpyrazine |
| | 3.653 | 2-ethyl-5-methylpyrazine |
| | 3.721 | 2-ethyl-6-methylpyrazine |
| | 3.72 | 2,3,5-trimethylpyrazine |
| | 3.91 | 2-ethyl-3,(5/6)-dimethylpyrazine |
| | 4.016 | 3,5-diethyl-2-methylpyrazine |
| | 4.654 | 2-methyl-(6/5)-acetylpyrazine |
| | 4.474 | 2-acetylpyrazine |
| **Carbonyls (C)** | 3.301 | 3-hydroxy-2-butanone |
| | 3.361 | 1-hydroxy-2-propanone |
| | 5.036 | 2-hydroxy-3-methyl-2-cyclopenten-1-one |
| | 5.219 | 3-ethyl-2-hydroxy-2-cyclopenten-1-one |
| **Furans and other oxygen-containing heterocyclic compounds (FU)** | 3.224 | dihydro-2-methyl-3(2H)-furanone |
| | 3.719 | 4-hydroxymethyl-2-methyl-1,3-dioxolane |
| | 4.681 | 5-methyl-2-furfurylalcohol |
| | 5.142 | 2,3-dihydro-5-hydroxy-6-methyl-4H-pyran-4-one |
| | 5.564 | **furaneol** |
| | 6.194 | 2,3-dihydro-3,5-dihydroxy-(2/6)-methyl-4H-pyran-4-one |
| | 6.756 | **5-hydroxymethy-2-furfural** [sic] |
| **Pyrroles and other N-containing heterocyclic compounds (PY)** | 5.381 | 2,3-dihydro-1-methyl-4(1H)-pyridinone |
| | 5.435 | 2-acetylpyrrole |
| | 5.879 | dimethyl-1H-benzimidazole |

`[D]` 10 PZ, 4 C, 7 FU, 3 PY = 24 identified compounds. **The "carbonyls" class contains no
α-dicarbonyl at all** — it is hydroxyketones and cyclotene-type enolones: > (p.2540) *"we have
supposed that hydroxyketones are formed in a similar way to α-dicarbonyls."* Methylglyoxal,
glyoxal, diacetyl are **[NEG] not listed**. The 1,3-dioxolane in FU is a glycerol–acetaldehyde
acetal (`[D]`, matrix artefact). HMF, furaneol and the two pyranones are lumped into FU with no
separate kinetics.

> (p.2540) "We have multiplied the total concentration for each class by a factor 2, to account
> for the fact that all volatiles have not been identified. This factor is, of course, completely
> arbitrary, and could be somewhat different."

---

## §3. EXPERIMENTAL SYSTEM — verbatim (p.2539) `[M]`

> "Twenty-four identical samples were prepared as follows: In a 50-mL tube (Kimble) with screw
> cap, **1.4 mmol of glucose and 1.45 mmol of alanine were dissolved in 10 g of phosphate-buffered
> glycerol (pH = 6.5)**. This mixture was heated at a given temperature for a given time. The
> reaction mixture was cooled to room temperature and dichloromethane with internal standards
> was added. The mixture was vortexed for 5 min and the dichloromethane layer was filtered over
> sodium sulfate into a GC vial. The samples were heated at **95, 110, 125, or 140 °C, for a total
> time of 180, 120, 60, and 30 min**, and taken out of the oven at time intervals of **30, 20, 10,
> and 5 min**, respectively."

| item | value | tag |
|---|---|---|
| sugar / amino acid | glucose 1.4 mmol; alanine 1.45 mmol (ratio 1 : 1.04) | [M] |
| matrix | 10 g phosphate-buffered **glycerol**, pH 6.5; buffer concentration **[NEG] not stated**; water content **[NEG] not stated** | [M] |
| nominal molarity | ≈ 0.18 M glucose, ≈ 0.18 M alanine (10 g glycerol ≈ 7.9 mL at ρ = 1.26 g/mL; buffer salts ignored) | [Z] |
| temperatures / durations | 95 °C / 180 min; 110 / 120; 125 / 60; 140 / 30 | [M] |
| sampling | 6 time points per temperature (30-min, 20-min, 10-min, 5-min steps) → 6 × 4 = 24 samples, one tube per point, **no replicates** | [Z] |
| heating lag | Eq. 2: **T = T_max + (T₀ − T_max)·exp(−t/τ)**, τ from measured T(t): **240 s at 95 °C to 273 s at 140 °C**; T₀ = room temperature | [M] |
| extraction | dichloromethane + internal standards (identity/amount **[NEG] not stated**), vortex 5 min, Na₂SO₄ | [M] |
| GC-MS (Table 3) | Analytical Applications Brielle; HP 6890 GC; WCOT fused silica CPWAX-52CB, 10 m × 0.10 mm, df = 0.2 µm; oven 40 °C (1 min) – 40 °C/min – 250 °C (9 min); HP 5973 MS; EI+, 70 eV; scan 30–300; sampling speed 2 → 9.4 scan/s | [M] |
| calibration | **[NEG] none described** — no response factors, no standards per compound; Fig. 10 y-axis is *"Relative concentration"* (unitless); Fig. 11 y-axis is *"Concentration (10⁻⁷ mol/L)"* | [M]/[NEG] |
| observed magnitudes (Fig. 11, raster) | PZ up to ≈ 5 × 10⁻⁷ M; FU ≈ 8 × 10⁻⁷ M; C ≈ 1.5 × 10⁻⁷ M; PY ≈ 0.9 × 10⁻⁷ M (140 °C, ×2 already applied) | [D] |
| implied yield | ≈ 10⁻⁶–10⁻⁵ of the 0.18 M glucose — extremely low; consistent with the "arbitrary" ×2 not being a calibration | [Z] |
| pyrazine detail (Fig. 10, p.2540) | methylpyrazine, 2-ethyl-6-methylpyrazine, 2-ethyl-3,(6/5)-dimethylpyrazine share one kinetic shape at all four T; **2-methyl-(6/5)-acetylpyrazine "follows markedly different kinetics"** (non-monotonic) — *"Our simplified kinetic modeling obviously cannot account for this difference, which points out its limits."* | [M] |
| solver | Eq. 1 *"solved numerically using Microsoft Excel, for all 4 temperatures simultaneously, using a temperature-dependent time-step"* | [M] |
| fitting | *"All rates were adjusted 'by hand,' starting from the average values determined in the previous section … only the most pleasing obtained after some 'playing around' with the rates."* | [M] |

---

## §4. EVERY RATE CONSTANT, PREFACTOR AND ACTIVATION ENERGY PRINTED

### §4.1 Table 1 (p.2535) — literature regressions, verbatim `[F]`

Title: *"Arrhenius rates for the different processes presented in Figure 1, determined from a
linear regression of lnR against 1/T"*. Footnote: *"ΔR₀ and ΔE correspond to the 95% confidence
interval on the regression. First-order processes are measured in per second; second-order are
measured in per second per molarity."*

| Process | Order | R₀ | ΔR₀ | E (kJ/mol) | ΔE |
|---|---|---|---|---|---|
| ARP formationᵃ | 2nd | 7.4 × 10¹⁴ | 2.1 × 10⁴ | 128.8 | 29.4 |
| AA loss | 1st | 2.6 × 10⁷ | 7.2 × 10⁰ | 84.1 | 5.8 |
| Sugar lossᵇ | 1st | 4.2 × 10⁶ | 3.3 × 10¹ | 75.7 | 10.9 |
| ARP decomposition | 1st | 2.3 × 10³ | 9.3 × 10¹ | 52.9 | 14.0 |
| HMF formationᶜ | 1st | 1.0 × 10¹¹ | 6.4 × 10⁴ | 109.3 | 31.3 |

ᵃ *"See Figure 2; data of Davies and others (1997) are excluded from the fit."* ᵇ *"See Figure 4;
data of Song and others (1966) and Bell and others (1998a) are excluded from the fit."* ᶜ *"See
Figure 6; data of Yaylayan and Forage (1991) are excluded from the fit."*

`[D]` ΔR₀ is a **multiplicative** factor on R₀ (ARP formation R₀ known to ×/÷ 2.1 × 10⁴),
inferred from the Fig. 8 min/max sets (§4.4). Sources per figure `[M]`: Fig. 2 ARP formation,
8 sources (Huyghues & Yaylayan 1995/1996, Lee 1984, Ge & Lee 1997, Baisier & Labuza 1992,
Warmbier 1976, Davies 1997, Leong & Wedzicha 2000); Fig. 3 ARP loss, 5 sources; Fig. 4 sugar
loss, 7 sources (*"scatter on 2 lines … the upper points gives an Ea of approximately 76
kJ/mol"*); Fig. 5 AA loss, 9 sources (*"from green peas … to cocoa beans … to model aqueous
systems"*); Fig. 6 HMF/FU, 10 sources, first-order only.

### §4.2 Table 2 (p.2537) — the "no fit" model set, verbatim `[F]`

Title: *"Arrhenius rates used to model Maillard chemistry."*

| Rate | Prefactor | E (kJ/mol) | R (20 °C) | R (100 °C) | R (180 °C) |
|---|---|---|---|---|---|
| R₁ /s | 1.7 × 10¹⁴ | 150.0 | 3.6 × 10⁻¹³ | 1.9 × 10⁻⁷ | 9.4 × 10⁻⁴ |
| R₂ /(Ms) | 7.4 × 10¹³ | 128.8 | 7.9 × 10⁻¹⁰ | 6.7 × 10⁻⁵ | 1.0 × 10⁻¹ |
| R₃ /s | 1.6 × 10⁶ | 73.3 | 1.4 × 10⁻⁷ | 8.7 × 10⁻⁵ | 5.7 × 10⁻³ |
| R₄ /s | 5.3 × 10³ (see C2) | 52.9 | 8.3 × 10⁻⁷ | 8.8 × 10⁻⁵ | 1.8 × 10⁻³ |
| R₅ /s | 1.0 × 10¹¹ | 109.3 | 3.1 × 10⁻⁹ | 4.8 × 10⁻⁵ | 2.4 × 10⁻² |
| R₆ /s | 2.2 × 10³ | 56.5 | 1.9 × 10⁻⁷ | 2.7 × 10⁻⁵ | 6.7 × 10⁻⁴ |
| R₈ /(Ms) | 3.6 × 10¹⁴ | 99.7 | 5.9 × 10⁻⁴ | 3.8 × 10⁰ | 1.1 × 10³ |
| R₉ /s | 1.0 × 10¹⁰ | 115.0 | 3.1 × 10⁻¹¹ | 7.7 × 10⁻⁷ | 5.4 × 10⁻⁴ |

`[Z]` All 24 R(T) cells recompute from (prefactor, E) with R = 8.314 J mol⁻¹ K⁻¹ to within
rounding **except** that R₄ requires 2.3 × 10³ (C2). R₇ and R₁₀ are absent by construction:
> (p.2538) "(1) The formation of FU from condensation of C by R₇ is supposed to be negligible,
> so that R₇ is always set to 0. (2) Strecker degradation is supposed to be a 'fast' process, so
> that R₉ and R₁₀ do not play any role."

Origin of each Table 2 entry `[M]`, from pp.2536–2538:

| rate | how the authors obtained it |
|---|---|
| R₁ | **no source given** (C5) |
| R₂ | Table 1 ARP formation (but 10¹³ vs 10¹⁴ — C1) |
| R₃, R₄ | Table 1 ARP decomposition (Ea 52.9, *"Ea of 55 kJ/mol"* in text) split by branch ratio: Huyghues-Despointes & Yaylayan 1995, **100 °C: amine-regenerating 0.35/h vs amine-losing 0.32/h**; Stahl & Parliament 1994, **≈ 4 at 200 °C** (amine-loss > regeneration). R₃'s 73.3 kJ/mol is *"estimated … from this analysis"* — no direct pyrrole data |
| R₅ | Table 1 HMF formation, as printed |
| R₆ | **Jusino, Ho & Tong 1997**, pyrazines in a **solid starch–lysine–glucose system, 80–120 °C, first-order, Ea 56.5 kJ/mol** — *"we use their data directly as an estimate of carbonyl formation"* (rate-determining-step argument) |
| R₈ | *"we use an Ea corresponding to the one measured for AA loss"* (84.1 → printed 99.7; see below) *"with a prefactor adjusted to reproduce the curve presented by Pripis-Nicolau and others (2000)"* (cysteine + methylglyoxal fully consumed **< 90 min at 25 °C**) |
| R₉ | Cremer & Eichner 2000a,b, Strecker aldehydes 80–110 °C, *"global Ea between 115 and 124 kJ/mol. We have used these Ea for R₉"* |

`[D]` **R₈'s printed Ea (99.7) does not equal Table 1's AA-loss Ea (84.1)** despite the sentence
quoted; 99.7 is also exactly the Ea later given to R₁₁. Unexplained. Other literature Ea quoted
in prose but not used `[C]`: Lerici 1990 ARP formation **160 kJ/mol**; Huang, Fu & Ho 1995
tetramethylpyrazine from 3-hydroxy-2-butanone 25–55 °C **Ea 79 kJ/mol**; Lerici 1990 CO₂ release
glycine/glucose 70–90 °C **102–115 kJ/mol**; Introduction survey of 66 sources: browning, HMF and
pyrazine Ea *"spread from 30 to 200 kJ/mol."*

### §4.3 Table 5 (p.2541) — "No fit" vs "Fit", verbatim `[F]`

Title: *"Rates used to fit the experimental concentration of PY, FU, SA, and PZ as a function of
time and temperature"* (SA is a title error, C3). Footnote: *"'Fast' indicates that the process
goes faster than the other steps, so that it does not contribute to the overall kinetics. The fit
is presented in Figure 11."*

| Process | No fit R₀ | No fit E (kJ/mol) | **Fit R₀** | **Fit E (kJ/mol)** |
|---|---|---|---|---|
| R₁ (/s) | 1.7 × 10¹⁴ | 150.0 | **0** | **0** |
| R₂ (/s/M) | 7.4 × 10¹³ | 128.8 | **5.0 × 10¹²** | **120.5** |
| R₃ (/s) | 1.6 × 10⁶ | 73.3 | **6.0 × 10¹** | **35.1** |
| R₄ (/s) | 2.3 × 10³ | 52.9 | **1.5 × 10⁵** | **52.9** |
| R₅ (/s) | 1.0 × 10¹¹ | 109.3 | **2.0 × 10¹¹** | **109.3** |
| R₆ (/s) | 2.2 × 10³ | 56.5 | **5.0 × 10⁵** | **66.5** |
| R₈ (/s/M) | 3.6 × 10¹⁴ | 99.7 | **5.0 × 10¹¹** | **83.1** |
| R₉ (/s) | 1.0 × 10¹⁰ | 115.0 | **1.0 × 10¹⁵** | **116.3** |
| R₁₀ (/s) | Fast | — | Fast | — |
| R₁₁ (/s) | — | — | **1.0 × 10¹⁰** | **99.7** |

(R₁₀ is printed with unit "/s" in Table 5 although Eq. 1 uses it bimolecularly — a units slip.)

`[Z]` Fit-set rate constants at **100 °C** (373.15 K), for prior-setting: R₂ 6.8 × 10⁻⁵ M⁻¹s⁻¹;
R₃ 7.3 × 10⁻⁴ s⁻¹; R₄ 5.9 × 10⁻³ s⁻¹; R₅ 1.0 × 10⁻⁴ s⁻¹; R₆ 2.4 × 10⁻⁴ s⁻¹; R₈ 1.2 M⁻¹s⁻¹;
R₉ 5.2 × 10⁻² s⁻¹; R₁₁ 1.1 × 10⁻⁴ s⁻¹. At 140 °C: R₄ 3.1 × 10⁻² s⁻¹, R₅ 3.0 × 10⁻³, R₆ 2.0 × 10⁻³,
R₁₁ 2.5 × 10⁻³ s⁻¹. The authors' own summary of the changes (p.2540–2541): *"R₃ (PY formation from
ARP) has smaller Ea; R₆ (C formation from RS) has larger Ea; the prefactor for R₄ (formation of RS
from ARP) is almost 2 orders of magnitude larger. These changes, however, remain well within the
limits of statistical uncertainty … The main difference, however, is to be found in the R₁ and R₁₁
terms … The plateauing observed in our model experiment could be fit only by the R₁₁ degradation
curve, which made the R₁ term unnecessary."*

### §4.4 Figure 8 caption — the CI-bounded alternative sets (p.2538) `[F]`

> "R₂ minimum, R₂ = 3.5 × 10¹⁰ exp(–99.4/RT); R₂ maximum, R₂ = 1.6 × 10¹⁹ exp(–158.2/RT);
> R₄ minimum, R₄ = 2.5 × 10¹ exp(–38.9/RT); R₄ maximum, R₄ = 2.2 × 10⁵ exp(–67.0/RT);
> R₅ minimum, R₅ = 1.6 × 10⁶ exp(–78.0/RT); R₅ maximum, R₅ = 6.5 × 10¹⁵ exp(–140.6/RT)"
> (activation energies in kJ/mol)

`[Z]` These are Table 1's E ± ΔE (128.8 ± 29.4 → 99.4/158.2; 52.9 ± 14.0 → 38.9/66.9; 109.3 ±
31.3 → 78.0/140.6) with prefactors slid along the regression (i.e. the CI is on the *slope*, and
the compensating intercept keeps the line pinned near the data centroid). The R₂ minimum/maximum
prefactors bracket 7.4 × 10¹⁴, not 10¹³ — further support for C1.

---

## §5. GOODNESS OF FIT — everything the paper says `[M]` / `[NEG]`

| claim | verbatim | quantitative content |
|---|---|---|
| Literature regressions | *"In all cases, the confidence interval is very large."* (p.2538) | Table 1 ΔE: ±5.8 to ±31.3 kJ/mol; ΔR₀ up to ×6.4 × 10⁴ |
| Fig. 7 prediction vs observation | *"There is good agreement, proving that the rates have been estimated well. Of course, we should not expect perfect agreement"* | log–log axes **10⁻¹⁴ to 10⁻⁴ M s⁻¹**, both axes; **[NEG] no R², no RMS log-deviation, no point count**. `[D]` from raster: HMF points lie within ≈ ±1 decade of the diagonal over 10⁻¹³–10⁻⁶; pyrazine points (×) sit systematically **above** the diagonal at observed 10⁻¹⁰–10⁻⁸ (over-predicted by ≈ 1–2 decades) |
| Sensitivity (Fig. 8) | *"The variations reach about 1 order of magnitude, which is not too much compared with the 9 orders of magnitude bandwidth of the data. We therefore conclude that small variations of the rate will not qualitatively change the results."* | ±1 decade on prediction from Table 1 CIs |
| Fig. 11 glycerol fit | *"Considering the high experimental error bars, the fit is rather good. It reproduces (1) the absolute magnitude of the volatiles; (2) the initial rate of volatile generation; and (3) the plateauing of the curves at long times."* | **[NEG] no R², no residuals, no error bars drawn, no replicates.** `[D]` raster: the 140 °C PZ and FU curves overshoot the last points; 110 °C FU is under-predicted by ≈ 2× at 60–120 min |
| Cross-checks claimed | *"the general browning term reproduces the order of magnitude of the observed browning rate in several model systems (Carabasa-Giribet and Ibarz-Ribas 2000), as well as the kinetics of formation of Strecker aldehydes (Cremer and Eichner 2000a, 2000b)."* | **[NEG] not shown** |
| Authors' own disclaimer | *"we do not claim that this is the best possible fit in terms of statistical accuracy … Clearly, it could be improved."* | — |

---

## §6. THE PSEUDO-ZERO-ORDER CORRELATION (FU and PZ vs initial T and concentration) `[M]`

> (p.2538) "we used the kinetic rates above to simulate the formation of furans and related
> compounds (FU) and of pyrazines (PZ), starting from the given temperature and concentration of
> reactants as given in several studies published in the literature. These experimental studies
> give the pseudo-zero-order rate of formation of the corresponding volatiles, which is the
> initial slope of the curve of flavor generation as a function of time. From our simulation, we
> measured this predicted initial slope and compare it with the experimental results in Figure 7."

Observed-rate sources (Fig. 7 caption): **HMF** — Lozano 1991; Mistry et al. 1995; Morales et al.
1995; Albalá-Hurtado et al. 1998; Gögüs et al. 1998; Bozkurt et al. 1999; Carabasa-Giribet &
Ibarz-Ribas 2000. **Pyrazines** — Huang, Bruecker & Ho 1989; Leahy & Reineccius 1989. Units:
**M s⁻¹** (initial slope). **[NEG] the paper prints no table of these observed rates, no T or
concentration per point, and no fitted correlation coefficient** — the "correlation" is the
Fig. 7 scatter only. Predicted rates use the Table 2 set (with R₂ at 10¹³, C1). Abstract wording:
*"the scheme was able to correlate the pseudo-zero-order rate of generation of FU and PZ (from the
literature) to the initial temperature and concentration of reactants."*

Additional stated conclusion (p.2541): *"the generation of pyrazines is optimized, as compared
with other flavor molecules and browning, when heating a glucose-alanine in a glycerol system at
140 °C for 30 min."* And on pH (not modelled): *"low pH favors the generation of furans and
derivatives, while high pH favors pyrazines. This could be easily accounted for, by making R₅
and/or R₆ pH-dependent."*

---

## §7. WHAT THE REPOSITORY CAN USE

**Topology mapping onto the trunk lane** (glucose/fructose/glycine → Schiff base → Amadori →
1-/3-deoxyosones → methylglyoxal, formic/acetic acid, HMF, DMHF, melanoidins):

| Jousse step | trunk analogue | what transfers | verdict |
|---|---|---|---|
| R₂ S + AA → ARP (2nd order) | glucose + glycine → Schiff base → Amadori, collapsed into one bimolecular step (cf. the repo's `amadori_over_schiff_pseudo_first_order`) | **Ea prior 128.8 ± 29.4 kJ/mol** (literature regression, 8 sources, aqueous/IM systems 37–110 °C, mostly glucose–lysine/glycine); fit value 120.5 in glycerol. k(100 °C) ≈ 6.7 × 10⁻⁵ M⁻¹s⁻¹ (model) or ≈ 6.9 × 10⁻⁴ (Table 1 line, C1) | **USE-Q as Ea prior; RATIO-ONLY for k** (10× ambiguity) |
| R₄ ARP → RS + AA | Amadori → 1-DG / 3-DG **with glycine release** | **Ea prior 52.9 ± 14.0 kJ/mol** (5 ARP-loss sources, 37–140 °C, kept unchanged by the fit). Branch: amine-regenerating ≈ amine-losing at 100 °C (0.35/h vs 0.32/h, morpholine-ARP) | **USE-Q** — the repo's Amadori → deoxyosone step should regenerate the amine; the 1-DG/3-DG split is **not** resolved here |
| R₅ RS → FU | 3-DG → HMF (and 1-DG → DMHF, since furaneol is in FU) | **Ea prior 109.3 ± 31.3 kJ/mol** (10 HMF sources, first-order, juices/milk/model, 20–120 °C); k(100 °C) 4.8 × 10⁻⁵ s⁻¹ (lit), 1.0 × 10⁻⁴ (fit) | **USE-Q as Ea prior for the 3-DG → HMF step only**; DMHF has no separate signal |
| R₆ RS → 2 C | deoxyosone retro-aldol → methylglyoxal (+ C₃ partner) | **Ea 56.5** is Jusino's *pyrazine* Ea in a *solid* starch system, relabelled; fit 66.5. The ×2 stoichiometry is the one structural fact | **PRIOR-ONLY, weak**; the stoichiometry is USE |
| R₁₁ volatiles → melanoidins | trunk melanoidin sink | single k for all classes: **1.0 × 10¹⁰ s⁻¹, 99.7 kJ/mol → 1.1 × 10⁻⁴ s⁻¹ at 100 °C** `[Z]`; needed to produce the plateau | **STRUCTURAL** (a uniform first-order loss of every volatile pool reproduces plateaus); the number is glycerol-specific |
| R₁ sugar direct degradation | caramelisation branch (glucose → 3-DG without amine) | 150 kJ/mol with no source; set to 0 by the fit | **REFUSE** |
| R₃ ARP → PY; R₈–R₁₀ Strecker/pyrazine | no trunk analogue | R₈ Ea 99.7 → 83.1; R₉ 115–124 (Cremer); Huang 1995 79 kJ/mol for pyrazine condensation | PRIOR-ONLY for a future pyrazine limb |
| Table 1 AA loss 84.1 ± 5.8; sugar loss 75.7 ± 10.9 | net glycine / glucose disappearance | the tightest CIs in the paper; span ~150 °C and wildly different matrices | **USE-Q as apparent-Ea sanity bounds** on the trunk's net reactant loss, not as step parameters |

**Absent from this paper for the trunk `[NEG]`:** formic acid, acetic acid, methylglyoxal,
glyoxal, diacetyl, any α-dicarbonyl; HMF- or DMHF-specific time series (both lumped into FU);
any pH or a_w dependence; any aqueous experiment of the authors' own; fructose; Schiff-base
reversibility.

**What a lumped "furans / carbonyls / pyrazines" surrogate would need from here:**
(i) the class membership lists of Table 4 (§2.4) as the definition of FU, C, PZ, PY; (ii) the
Eq. 1 topology with R₇ = 0, R₁₀ fast (so PZ ≈ ½ × flux through R₈ net of R₉), the ×2 on R₆, the
amine-regeneration on R₄, and the uniform R₁₁ sink; (iii) the Table 5 "Fit" (R₀, Ea) set as an
*initialisation*, with the 100 °C values in §4.3 `[Z]`; (iv) the ×2 unidentified-compound
multiplier made explicit as a free scale; (v) the Eq. 2 heating lag (τ 240–273 s) if the 95–140 °C
Fig. 11 curves are ever digitised as a hold-out — they would need to be read from the raster
(no table exists) and are in 10⁻⁷ mol/L after ×2.

## §8. CAVEATS

1. **Glycerol, not water.** The only experiment is in phosphate-buffered glycerol (a_w ≈ 0, water
   content unstated). The authors: *"The rates given in Table 5 are valid, strictly speaking, only
   for a model glucose-alanine system in glycerol at a pH of 6.5"*; *"it is well known that the
   water activity affects the kinetics in a nontrivial way … The nature of the solvent and its
   physical state also play important roles."* "pH 6.5" in glycerol is an operational, not a
   thermodynamic, quantity. The 1,3-dioxolane in FU is a glycerol artefact.
2. **Lumped classes.** PY/FU/C/PZ sum unequal compounds; "carbonyls" are hydroxyketones and
   cyclic enolones, not dicarbonyls; HMF, furaneol and pyranones share one rate. The authors'
   central postulate — *"the kinetics of formation for all compounds in the same class is
   similar"* — is contradicted by their own acetylpyrazine data (Fig. 10).
3. **Hand fit, no statistics.** Excel, "by hand", "most pleasing", ×2 arbitrary scaling, single
   tubes, no calibration described, no R²/residuals/CIs on Table 5. Prefactors are
   compensating quantities (R₄ moved ×65, R₉ ×10⁵ between "no fit" and "fit").
4. **Literature regressions are cross-system.** 66 sources; peas, cocoa, milk, juices, powders,
   aqueous models; points excluded when "clearly off-trend" (the authors admit *"we make a choice
   and therefore bias the resulting rates"*); CIs of ±14–31 kJ/mol on the three step Ea. Treat as
   priors with those widths, never as measured trunk Ea.
5. **Internal number defects (C1, C2).** R₂ prefactor 10¹⁴ vs 10¹³ and R₄ 5.3 vs 2.3 × 10³ are
   printed inconsistently; §0.2 resolves which value each table actually used, but the paper
   never does.
6. **Borrowed and relabelled Ea.** R₆ (carbonyl formation) = Jusino's pyrazine Ea in a solid
   system; R₈ = "AA loss" Ea but printed as 99.7 ≠ 84.1; R₉ = Cremer's Strecker-aldehyde Ea;
   R₁ has no source; R₃ Ea is inferred from a branch ratio at 100 °C and 200 °C.
7. **Verification is order-of-magnitude only.** Fig. 7 spans nine decades; the authors' own
   sensitivity is ±1 decade; pyrazines are systematically over-predicted (§5, `[D]`).
8. **Scope of the model itself.** No reversibility, no acids, no water, no pH, no a_w, no
   Heyns/ARP distinction, no 1-DG/3-DG split, single melanoidin sink with one k for five pools.

# Szajewski & Whitesides 1980 — EXTRACTION (rates of GSSG reduction by 9 mono-/dithiols in 0.066 M phosphate pH 7.0, 30.0 °C, argon, glyoxalase-I-coupled assay; a three-parameter Brønsted model for thiolate–disulfide interchange fitted to ~50 literature rates; equilibrium constants for dithiol + oxidised lipoamide / GSSG; separation of literature K_obsd into thiolate (K^S−) and thiol (K^SH) parts)

### THE STANDARD LINEAR FREE-ENERGY MODEL FOR THIOL–DISULFIDE EXCHANGE: **log k (M⁻¹ min⁻¹) = 7.0 + 0.50 pKa(nuc) − 0.27 pKa(central S) − 0.73 pKa(leaving group)** (eq 36; r² 0.98 over ~50 rates; aryl-aryl: C = 7.9), with two equally good alternatives (eqs 37-38). For GSSG specifically: **log k1(RS⁻ + GSSG) = −1.29 + 0.50 pKa(RSH)** (r² 0.89; pH 7.0, 30 °C). Equilibria: thiolate exchange **log K^S⁻ ≈ 1.21 × 2(pKa^RSH − pKa^ESH)** (double exchange; β_nuc − β_lg ≈ 1.2, ±20 %), thiol exchange K^SH ≈ 1 with steric scatter. **One temperature (30 °C), no Ea, ionic strength not controlled.**

**Source on disk:** `data/articles/szajewski1980.pdf` (16 pp., ACS scan with OCR text layer; noisy). Equations and
tables taken from `pdftotext -layout` and checked against the context of the surrounding text; I did not render the
figures (Figs. 3, 5, 7, 10 are plots whose fitted lines are printed as equations).

## 0. Identity

| field | value |
|---|---|
| Title | "Rate Constants and Equilibrium Constants for Thiol-Disulfide Interchange Reactions Involving Oxidized Glutathione" |
| Authors | Richard P. Szajewski, **George M. Whitesides** — MIT |
| Venue | J. Am. Chem. Soc. **102** (6), 2011-2026, 12 Mar 1980. Received 8 May 1978 |
| DOI | not printed on the 1980 page; the brief's **10.1021/ja00526a042** matches title, journal, volume, pages. **Identity verified.** |

## 1. Why it matters

SULFUR_FLUX_MAP sec. 5b: in Zhou 2023c (Ala-xylose ARP + glutathione, pH 7.5, 120 °C) the measured MFT collapses while the
model's keeps rising; ~81 % of the cysteine moiety is unaccounted for, probably GSSG and mixed disulfides, which are
thiol–disulfide exchange partners for MFT that the engine lacks. A rate law for RS⁻ + R'SSR'' with pKa dependence is what
such a step needs.

## 2. Methods as they matter (pp. 2011-2014, 2023-2024)

- **Rates.** GSSG 0.054-0.64 mM (standard 0.32-0.35), thiol 0.49-14.5 mM (standard ~5 mM), methylglyoxal 0.17-1.50 mM,
  glyoxalase-I 0.3-24 U/mL, **0.066 M phosphate pH 7.0, 30.0 ± 0.5 °C, argon**. GSH released is trapped as S-lactoyl-GSH
  (240 nm, ε 3.37 mM⁻¹ cm⁻¹). Integrated second-order law (eq 13); 10-60 % conversion. Corrected for a slow
  phosphate-catalysed methylglyoxal + thiol blank. Reproducibility ±15 %; typically triplicate.
- **Thiolate rate constant** from the observed one: k1 = k1_obsd (1 + 10^(pKa^RSH − pH)) (eq 10).
- **Equilibria.** Dithiol + oxidised lipoamide with lipoamide dehydrogenase/NAD⁺, 0.05-0.2 M phosphate pH 7.0, 30 °C,
  argon; converted to K vs GSSG and to E°' (vs SHE, pH 7, 30 °C).
- **Ionic strength:** not stated / not controlled (phosphate only). **Temperature:** 30 °C only.

## 3. Numbers

### Table III (p. 2014), RSH + GSSG, pH 7.0, 30 °C (M⁻¹ min⁻¹; ±15 %)

| thiol | pKa (2nd pKa) | k_obsd | k1 (thiolate) |
|---|---|---|---|
| glycol dimercaptoacetate (GMA) | 7.7 (9.0) | 90 | 4.9e2 |
| 1,3-dimercapto-2-propanol (DMP) | 9.0 (10.3) | 9.3 | 9.3e2 |
| dithiothreitol | 9.2 (10.1) | 14.1 | 2.2e3 |
| DTT/DTE mixture | 9.2 (10.1) | 8.4 | 1.1e3 |
| 1-thioglycerol | 9.56 | 7.9 | 2.5e3 |
| N-acetylcysteine | 9.5 | 5.2 | 1.8e3 |
| 2-mercaptoethanol | 9.6 | 8.7 | 3.4e3 |
| thioglycolic acid | 9.8 | 7.3 | 4.6e3 |
| 3-mercaptopropionic acid | 10.6 | 3.2 | 1.2e4 |

Fig. 3 lines: **log k1 = −1.29 + 0.50 pKa** (β_nuc = 0.50, r² 0.89); log k1_obsd = −1.29 + 0.50 pKa − log(10^(pKa−7.0) + 1).
GSH pKa taken as 8.72. **Cysteine and cysteamine could not be measured** (they react rapidly with methylglyoxal in the
assay); thioacetic acid indistinguishable from blank. Tables I-II: k_obsd invariant (±~15 %) to 10× changes of enzyme,
thiol, methylglyoxal and GSSG (2-mercaptoethanol 7.3-9.5; DTT 8.4-14.6 M⁻¹ min⁻¹; DTT is lowered at 0.2× enzyme and 0.5×
methylglyoxal).

### General Brønsted models (p. 2021; k in M⁻¹ min⁻¹, thiolate rate constant; fitted to rates from this work and refs. 12, 13, 21, 22, Gorin 1967 — mostly 25-30 °C, aqueous, pH 6-8)

| eq | model | r² | constraint |
|---|---|---|---|
| 36 | log k = **7.0 + 0.50 pKa^nuc − 0.27 pKa^c − 0.73 pKa^lg** | 0.98 (slope 0.98, intercept 0.11 vs expt) | β_nuc fixed 0.5; β_c + β_lg = −1.0 |
| 37 | log k = 5.5 + 0.63 pKa^nuc − 0.30 pKa^c − 0.63 pKa^lg | 0.96 | β_c = −0.30, β_nuc = −β_lg |
| 38 | log k = 6.3 + 0.59 pKa^nuc − 0.40 pKa^c − 0.59 pKa^lg | 0.96 | β_nuc = −β_lg, β_c free |
| 39 | log k_obsd = log k − log(1 + 10^(pKa^nuc − pH)) | — | (printed without the "log"; misprint) |

If nucleophile and leaving group are both aromatic thiols, C = 7.9 in eq 36. Thioglycolic acid / methyl thioglycolate
are outliers. Authors: "the values of these Brønsted coefficients are not sharply defined"; eqs 36 and 38 are used
interchangeably as kinetic models.

### Equilibria

- **Thiolate exchange:** Fig. 7B line **log K^S⁻ = 2.42 (pKa^RSH − pKa^ESH)** through (0,0) for 2RS⁻ + ESSE ⇌ RSSR + 2ES⁻
  (alkyl thiols 1-8 of Table IV), i.e. slope 1.21 against 2ΔpKa → **β_nuc − β_lg ≈ 1.2 (±≥20 %)**. Per single exchange
  (RS⁻ + ESSE ⇌ RSSE + ES⁻) this is log K ≈ 1.2 ΔpKa (mine, assuming the two steps are equivalent).
- **Thiol exchange:** K^SH shows no pKa correlation; ~0.1-30, lowered by steric hindrance (Table IV).
- **Table IV (selected, dimensionless):** cysteamine + cystine pH 7.4: K_obsd 3.6; **GSH + cystine pH 6.0: K_obsd 3.0,
  K^S⁻ 28, K^SH 3.0**; cysteamine + GSSG pH 7.4: 1.7; thioglycolate + cystine pH 6.0: 3.3; **cysteine + an aryl disulfide
  (BSSB) pH 6.0: K_obsd 0.98, K^S⁻ 1.2e4, K^SH 0.12**; thiophenol + dialkyl disulfide (hydrocarbon) 0.16.
- **Table V:** GSH pKa 8.7, E°' −0.205 V; 2-mercaptoethanol −0.207 V; DTT −0.323 V (pH 7.0, 30 °C, vs SHE).
- Text: equilibria for reduction of **diaryl disulfides by aliphatic thiols "lie so far toward aryl thiol"** that they
  could not be measured.
- Side datum (eq 5): **GSH + methylglyoxal ⇌ hemithioacetal, K5 = 3 mM** (dissociation constant, cited ref. 14).

## 4. Numbers the repository can use

| step | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| RS⁻ + R'SSR'' → RSSR' + R''S⁻ | Brønsted model | eq 36 (or 38) | M⁻¹ min⁻¹ | aq., ~25-30 °C, pH 6-8 | p. 2021 | fitted LFER over ~50 literature rates |
| RS⁻ + GSSG | Brønsted | log k1 = −1.29 + 0.50 pKa | M⁻¹ min⁻¹ | 0.066 M phosphate pH 7.0, 30 °C | Fig. 3, Table III | measured (9 thiols) |
| thiolate exchange equilibrium | β_nuc − β_lg | ≈ 1.2 | — | aq. | Fig. 7B | fitted, ±20 % |
| thiol exchange equilibrium | K^SH | ~1 (0.1-30) | — | aq. | Table IV | literature compilation |
| GSH + cystine | K_obsd | 3.0 | — | pH 6.0 | Table IV row 6 | secondary (Eldjarn & Pihl) |
| GSH + methylglyoxal hemithioacetal | K_d | 3 | mM | assay conditions | eq 5 | secondary |
| any Ea / ΔH | — | **NOT REPORTED** | — | — | — | — |

## 5. What this can and cannot anchor

**Can:** a pKa-based rate and equilibrium model for thiolate–disulfide interchange among alkyl thiols at ~30 °C, including
GSSG; the pH dependence through the thiolate fraction (eq 39).
**Cannot:** temperature dependence (none); heteroaryl thiols such as MFT/FFT (no furan thiols; pKa's of MFT/FFT are not in
this paper); cysteine itself as nucleophile on GSSG (not measurable in the assay); ionic-strength effects.

## 6. Flags

1. Eq 39 misprint (missing "log").
2. The Brønsted coefficients are not uniquely determined (eqs 36-38 fit equally); only the qualitative charge picture and
   the predicted rates are robust (authors' own statement).
3. The eq 36 data set mixes conditions (25-30 °C, various buffers and pH).
4. The OCR text layer is noisy; Table V columns are partly garbled — only the values quoted above were checked against
   their textual context.

## 7. Use in the model

- **An MFT sink via glutathione disulfides (Zhou 2023c).** A step MFT⁻ + GSSG ⇌ MFT-SG + GS⁻ (and MFT⁻ + mixed disulfides)
  can be written with eq 36/38 for k and the pKa relation for K, both **declared** at 30 °C with a declared Ea band (none is
  measured here). Thermodynamic caution from this paper: for a thiol more acidic than GSH (aryl-like), the thiolate
  exchange equilibrium lies toward the free, more acidic thiolate (log K ≈ 1.2 (pKa_MFT − 8.7) per exchange; with an
  *assumed* MFT pKa of ~6.5-7, K ≈ 10⁻² to 10⁻³), and aryl thiols stay reduced against alkyl disulfides. So GSSG is a
  **fast but thermodynamically weak** MFT sink at equilibrium; a large MFT-SG pool needs GSSG ≫ MFT, or an irreversible
  follow-up. Whether that explains Zhou's 59 → 8 µg/L collapse is a modelling question; this paper alone does not settle
  it. MFT's pKa must come from elsewhere.
- **Cysteine / GSH / cystine exchange** (Table IV: GSH + cystine K 3.0 at pH 6) is available if the engine ever carries
  cystine or GSSG explicitly.
- No bearing on the furfural barrier swap or B61.

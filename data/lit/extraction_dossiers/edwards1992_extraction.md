# Edwards & Wedzicha 1992 — EXTRACTION (3-deoxyhexosulose (DH) 10-50 mM + mercaptoethanol, glutathione, N-acetylcysteine, sulfite; and + cysteine / homocysteine / cysteamine 25 mM; 50 mM acetate pH 5.5, 0.5 mM EDTA, ± glycine 0-1.0 M, 55 °C; thiol by Ellman; a cabbage-slurry test)

### THE FOLLOW-UP TO WEDZICHA & EDWARDS 1991: it re-states the DH + thiol rate law (k1, k2) for four nucleophiles in one table, and its aminothiol figure carries the decisive control — CYSTEAMINE, which has no carboxyl and cannot undergo a Strecker decarboxylation, also binds DH fast and then slowly liberates thiol. So "fast thiazolidine, slow thiol release" is a property of aminothiols with DH, not proof of a Strecker step.

**Source on disk:** `data/articles/edwards1992.pdf` (10 pp.: a Taylor & Francis cover + journal pp. 461-469, typeset
scan). Read from 150-dpi renders through the Read tool; **Fig. 5 digitised from a 300-dpi render** (scratchpad
`img/edwF5b.png`), tagged **[D]**. Table 1 re-typed from the render. Figs. 2-4 (NAC stability, HPLC traces) and Fig. 6
(cabbage) were read only coarsely.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics and mechanism of the reaction between 3-deoxyhexosulose and thiols" ✔ |
| Authors | A. S. Edwards, B. L. Wedzicha — Procter Department of Food Science, University of Leeds ✔ |
| Venue | Food Additives and Contaminants **9**(5), 461-469 (1992) ✔ |
| DOI as printed | **10.1080/02652039209374098** ✔ matches the brief |
| Companion | `wedzicha1991_extraction.md` (the full kinetic paper: rate law derivation, pH scan, D2O, NaBH3CN, cysteine binding at 10 °C, TLC of the liberated thiol). This paper cites it for the cysteine 1:1 stoichiometry and the Strecker proposal. |

## 1. Why it matters

Same gaps as the companion — (b) cysteine + dicarbonyl ⇌ thiazolidine, slow release; (a) whether the slow step is a
Strecker giving mercaptoacetaldehyde (MAA). What this paper adds:

1. **One table of k1 / k2 / stoichiometry for four nucleophiles** at 55 °C, pH 5.5 (Table 1), including NAC and sulfite.
2. **The aminothiol comparison** (Fig. 5, no glycine): cysteine, homocysteine and cysteamine all bind DH within ~0.5-1.5 h at
   55 °C and then liberate thiol over 2-7 h. **Cysteamine cannot decarboxylate**, yet liberates ~3.3 mM; cysteine
   liberates ~5.8 mM over the same window. The authors nonetheless assign the slow step for cysteine and homocysteine to
   a Strecker degradation (p. 468) and leave cysteamine's release unexplained.
3. **Stoichiometries with pyruvaldehyde** (methylglyoxal, an engine species `MGO`): DH : thiol 1 : 1.01 (cysteine),
   1 : 1.02 (cysteamine), 1 : 0.86 (homocysteine); cysteamine + DH at 25 °C 1 : 0.94, "no subsequent reaction" at 25 °C.
   Pyruvaldehyde binds aminothiols "more rapidly" than DH (not quantified).
4. A **"'meaty' aroma quickly develops with cysteine"**, even at room temperature (p. 468).

## 2. Methods as they matter to a model

- **Charge.** "DH (10-50 mM), thiol (10-50 mM) and glycine (0-1·0 M)" in "pH 5·5 acetate buffer (50 mM) containing EDTA
  (0·5 mM)" (p. 462). 2 mL vials, sealed, **55 °C**. Thiol by Ellman's reagent (412 nm). Stoichiometry from residual thiol
  after completion with [thiol] > [DH].
- **Crotonaldehyde control.** 20 mM crotonaldehyde + 50 mM ME, 2 h, 55 °C, pH 5.5: thiol plateau 30.5 ± 1.2 mM (n = 6) without
  glycine, 31.3 ± 1.0 mM (n = 6) with 0.5 M glycine → **1:1** thiol addition to an α,β-unsaturated carbonyl.
- **NAC stability.** 50 mM NAC in acetate pH 5.5 at 55 °C, followed by Ellman ± Dowex 50W-X8 (H⁺) ion exchange and HPLC
  (Spherisorb 5ODS, 25 mM KH2PO4/Na2HPO4 pH 7.0, 240 nm).
- **Aminothiols (Fig. 5).** [aminothiol] = 25 mM, [DH] = 20 mM, 50 mM acetate pH 5.5, 55 °C, **no glycine listed**.
- **Cabbage.** Freeze-dried blanched cabbage 2 g in 25 mL of 10 mM thiol or sulfite at 55 °C; filtrate + washings by Ellman.
  Cabbage background ≈ 0.4 mM thiol-equivalent.
- **Not reported:** any temperature other than 55 °C (and 25 °C for one stoichiometry), any pH other than 5.5 for the
  aminothiols, any product identification, H2S, error bars on Fig. 5.

## 3. Tables re-typed and figures digitised

### Table 1 (p. 463). k1, k2 and stoichiometries for thiols and sulfite with DH at 55 °C and pH 5.5

| nucleophile | k1 / h⁻¹ | k2 / M⁻¹ | stoichiometry DH : thiol |
|---|---|---|---|
| mercaptoethanol | 0.029 | 1.31 | 1 : 2 |
| glutathione | 0.020 | 1.35 | 1 : 2 |
| N-acetylcysteine | 0.007 | 1.84 | 1 : < 1 |
| sulphite | 0.005 | 4.45 | 1 : 1 |

Rate law (p. 462): **initial rate = k1(1 + k2[glycine])[DH]**, first order in DH, zero order in thiol, glycine catalysis
first order. ME and GSH values are those of Wedzicha & Edwards 1991 Table 1; NAC is new here.

### Fig. 5 [D] (p. 467). Total thiol (mM) vs time (h); [aminothiol] = 25 mM, [DH] = 20 mM, acetate 50 mM pH 5.5, 55 °C

Read error ± 0.1 mM, ± 0.05 h. t = 0 (25 mM) is the charge, not plotted.

| t (h) | 0.22 | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 | 4.0 | 5.0 | 7.0 |
|---|---|---|---|---|---|---|---|---|---|
| cysteine (○) | 11.2 | 10.2 | 10.5 | 10.7 | 11.4 | 11.9 | 13.2 | 14.1 | 15.95 |
| homocysteine (□) | 10.5 | 7.8 | 7.9 | 8.0 | 8.0 | 8.4 | 9.1 | 9.6 | 9.5 |
| cysteamine (△) | 11.6 | 10.5 | 9.9 | 9.45 | 9.7 | 10.1 | 11.0 | 11.5 | 12.7 |

### Fig. 2 (p. 464), coarse [D] (± 0.3 mM): NAC 50 mM, acetate pH 5.5, 55 °C

Total thiol ≈ 48.8 mM at ~25 h → 43.2 mM at ~340 h; thiol not retained by the cation exchanger 46.6 → 38.8 mM. HPLC
(text): NAC **43 mM after 93 h** vs 42.5 mM non-retained thiol. In 0.1 M HCl at 55 °C NAC converts to cysteine (Fig. 4).
NAC slowly isomerises (proposed via S-acetylcysteine) at pH 5.5; only traces of cysteine form.

### Fig. 6 (p. 468): not digitised

10 mM thiol or sulfite in cabbage slurry, 55 °C, to 140 h. Cysteine is lost fastest; at ~140 h ME, GSH and
mercaptoethanesulfonic acid have reacted about **twice** as much as sulfite (text). Not used.

### Arithmetic (mine)

**1. Binding extent at 55 °C.** Minimum thiol: cysteine 10.2 mM at 0.5 h → **14.8 mM bound = 74 % of DH**;
homocysteine 7.8 mM at 0.5 h → 17.2 mM bound (86 % of DH); cysteamine 9.45 mM at 1.5 h → 15.55 mM bound (78 %), and
**cysteamine binds more slowly** (still falling at 1.5 h, where cysteine has already turned).

**2. Liberation.** From the minimum to 7 h:

| aminothiol | liberated (mM) | window (h) | mean rate (mM h⁻¹) | pseudo-first-order on bound (h⁻¹) | shape |
|---|---|---|---|---|---|
| cysteine | 5.75 | 0.5 → 7 | 0.88 | **0.076** | sigmoidal, still rising at 7 h |
| homocysteine | 1.7 | 2 → 5 | 0.57 (then flat) | ~0.035 (2-5 h) | plateau 9.5-9.6 from 5 h |
| cysteamine | 3.25 | 1.5 → 7 | 0.59 | **0.044** | sigmoidal, plateau ~12.7 |

Cysteine's 0.076 h⁻¹ **reproduces the companion's ~0.07 h⁻¹** (pH 5.5 with 0.5 M glycine, `wedzicha1991_extraction.md`
arithmetic 3) without glycine — glycine does not drive the liberation step.

**3. The cysteamine control.** Liberation of 58 % as much thiol as cysteine, at 58 % of the pseudo-first-order rate, from
an aminothiol that has no carboxyl. If cysteine's liberation were purely Strecker (decarboxylative), cysteamine should
show none. Two non-exclusive readings: (i) part of every aminothiol's release is a non-Strecker process (slow thiazolidine
ring-opening as DH converts to other products, or transfer to a less stable adduct), and the Strecker share of cysteine's
release is at most ~40 % (5.75 − 3.25 = 2.5 mM of 5.75, if the non-Strecker part is the same size as cysteamine's);
(ii) cysteamine undergoes a different, non-decarboxylative transamination-like cleavage. **The paper does not resolve it.**
Wedzicha 1991's TLC (liberated thiol ≠ cysteine) is not repeated for cysteamine.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | evidence class |
|---|---|---|---|---|---|
| DH activation (RDS), ME / GSH / NAC / sulfite | k1 | 0.029 / 0.020 / 0.007 / 0.005 | h⁻¹ | 55 °C, pH 5.5 acetate | **measured** (Table 1) |
| glycine catalysis | k2 | 1.31 / 1.35 / 1.84 / 4.45 | M⁻¹ | same | measured (Table 1) |
| α,β-unsaturated carbonyl + ME | stoichiometry | 1 : 1 | — | crotonaldehyde, 55 °C, 2 h | measured |
| aminothiol + pyruvaldehyde | stoichiometry | Cys 1:1.01; cysteamine 1:1.02; homoCys 1:0.86 | — | 55 °C (implied), pH 5.5 | measured (text) |
| cysteamine + DH, 25 °C | stoichiometry | 1 : 0.94, no further reaction | — | 25 °C, pH 5.5 | measured (text) |
| Cys + DH → adduct | extent at 0.5 h | 74 % of DH | — | Cys 25 / DH 20 mM, 55 °C, no gly | derived_from_digitised |
| Cys–DH adduct → thiol | pseudo-first-order | **0.076** | h⁻¹ | same | **derived_from_digitised** |
| cysteamine–DH adduct → thiol | pseudo-first-order | **0.044** | h⁻¹ | same | derived_from_digitised |
| homocysteine–DH adduct → thiol | pseudo-first-order | ~0.035 | h⁻¹ | same, 2-5 h only | derived_from_digitised |
| NAC at pH 5.5, 55 °C | loss | ~50 → 43 mM in 93 h | mM | acetate | measured (HPLC, text) |

No Ea, no pH series for aminothiols, no H2S, no product identification.

## 5. What this can and cannot anchor

**Can:**
- The ME/GSH/NAC/sulfite k1, k2 as **measured constants for a hexose 3-deoxyosone's thiol-independent activation** at
  55 °C (relevant if the engine ever carries a 3-deoxyosone + thiol (MFT/FFT) sink — `ch_add_mft_tdp`, `ch_add_fft_tdp`
  exist and use the generic `k_add`). The rate law says such a sink is **zero order in the thiol** and limited by the
  osone's own conversion, with 2 thiols consumed per osone.
- A second, glycine-free estimate of the cysteine–DH adduct's liberation step: **0.076 h⁻¹ at 55 °C, pH 5.5**.
- The **cysteamine control** that caps the Strecker share of that liberation.
- Pyruvaldehyde (MGO) binds cysteine 1:1 — the same thiazolidine chemistry for the engine's `MGO` (+ `r_cys_actz`).

**Cannot:**
- No temperature dependence; no pentose osone; no H2S; the liberated thiol's identity is not established here.
- The aminothiol figure has no error bars and only one charge ratio.

## 6. Flags

1. **The authors' Strecker assignment is under-determined by their own data** (cysteamine liberates thiol too; arithmetic 3).
2. Table 1's ME/GSH values duplicate the 1991 paper; do not double-count them as independent measurements.
3. Fig. 5 first points (0.22 h) sit above the minima — binding is not instantaneous at 55 °C: cysteine ~0.5 h, cysteamine
   ~1.5 h to the minimum.
4. NAC is unstable at pH 5.5 / 55 °C (~14 % thiol loss in ~340 h): irrelevant to the engine but a warning for NAC-fed controls.
5. Registry: paper absent from `data/keys/papers.yml`; DH, MAA, homocysteine, cysteamine are not engine species
   (cysteamine is relevant to Hofmann 1995b's 2-acetyl-2-thiazoline route, `hofmann1995b_extraction.md`).

## Use in the model

- **Gap (b).** With the companion: cysteine + 3-deoxyosone binds fast (≤ 0.5 h at 55 °C, ~74 % of the osone at 25:20 mM)
  and the adduct transforms at **~0.07-0.08 h⁻¹ (55 °C, pH 5.5), glycine-independent**. That is a slow, first-order
  conversion of a cysteine-osone adduct — the kind of step B59's de Roos path needs between "bound" and "released". Carry
  it as a **declared band with no Ea** (hexose analogue, 55 °C).
- **Gap (a).** **Weakens** the inference "slow release = Strecker = MAA": cysteamine, which cannot Strecker, releases
  thiol at ~60 % of cysteine's rate. At most ~40 % of cysteine's liberated thiol is attributable to a decarboxylative
  (Strecker) process on this evidence. Any MAA yield drawn from these two papers should be bounded accordingly.
- **Gap (f), cysteamine.** Cysteamine binds an α-dicarbonyl (DH, pyruvaldehyde) 1:1 like cysteine — any cysteamine pool
  in a cysteamine/thiazoline lane would be partly sequestered by dicarbonyls; relevant to Hofmann 1995b, not quantified here.
- **Caveats:** 55 °C only; hexose 3-deoxyosone; Ellman thiol is a sum over thiols; no product identification.

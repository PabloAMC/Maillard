# Wedzicha & Edwards 1991 — EXTRACTION (3-deoxyhexosulose (DH) 5-50 mM + mercaptoethanol / glutathione / cysteine 5-50 mM, ± glycine 0-1.0 M, 50 mM acetate pH 5.5 (pH 3.85-6.00 scan), 0.5 mM EDTA, 55 °C and 10 °C; one pH 6.9 phosphate arm ± NaBH3CN; thiol by Ellman / iodometry)

### THE FULL PAPER BEHIND EDWARDS 1992 AND THE ONLY KINETIC RECORD IN THE CORPUS OF CYSTEINE MEETING A 3-DEOXYOSONE: cysteine binds DH 1:1 within minutes (zero order in cysteine, ~3.8 × 10⁻⁴ M min⁻¹ at 24 mM DH even at 10 °C), and over hours at 55 °C the adduct liberates a NEW thiol that is not cysteine (TLC Rf 0.12 vs 0.27), at ~0.9-2 mM h⁻¹ — the authors propose mercaptoacetaldehyde from a Strecker degradation of the thiazolidine. For mercaptoethanol and glutathione the rate law is first order in DH, zero order in thiol and H⁺ (pH 3.85-6.00), with k = k1(1 + k2[glycine]).

**Source on disk:** `data/articles/wedzicha1991.pdf` (16 pp., the full article pp. 71-86; an image-only scan with an
OCR layer, Acrobat Capture 2003). Read from 150-dpi page renders through the Read tool; **Figs. 3 and 4 digitised
from 300-dpi renders** (scratchpad `img/wedF3.png`, `img/wedF4.png`). Every digitised value is tagged **[D]**.
Tables 1-3 are re-typed from the renders. No supplementary material.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics of the Reaction of 3-Deoxyhexosulose with Thiols" ✔ |
| Authors | B. L. Wedzicha, A. S. Edwards — Procter Department of Food Science, University of Leeds ✔ |
| Venue | Food Chemistry **40** (1991) 71-86; received 6 Dec 1989, revised and accepted 10 May 1990 ✔ |
| DOI | **not printed** on the article; the PDF's metadata title is "PII: 0308-8146(91)90020-O", which is the brief's DOI suffix 10.1016/0308-8146(91)90020-o ✔ |
| Companion | `edwards1992_extraction.md` (Edwards & Wedzicha 1992, Food Addit. Contam.: same lab, same protocol; adds N-acetylcysteine, cysteamine, homocysteine, pyruvaldehyde stoichiometries and a cabbage test) |
| Lineage | Wedzicha & Kaban 1986 (DH + S(IV) + glycine kinetics, Food Chem. 22:209), whose rate law this paper parallels |

## 1. Why it matters

The engine's 3-deoxyosone of a pentose (`TDP`) has **no reaction with cysteine** at all (`sulfur.py`: TDP is made by
`r_pent_tdp`/`r_pent_caramel`/ARP routes and leaves only to furfural, osone decay, melanoidin and thiol adducts).
Cysteine's only dicarbonyl partner is the 1-deoxyosone `DPO` in the Strecker step `r_dpo_ddp`, whose sulfur goes to
inert `FRAG_S`. This paper — on the hexose analogue DH (3-deoxyglucosone) — shows:

1. **Fast, 1:1 binding** of cysteine to the 3-deoxyosone: at 10 °C, pH 5.5, no glycine, the initial rate is
   **3.8 × 10⁻⁴ M min⁻¹ at [DH] = 24.4 mM and independent of [cysteine] over 10-50 mM** (Table 3). The product is
   proposed to be the thiazolidine (or, less likely, a thiazine); its thiol is not released by the Ellman assay.
2. **Slow liberation of a different thiol** at 55 °C (Fig. 3, Fig. 4): a DTNB-positive spot at Rf 0.12, not cysteine
   (Rf 0.27). The authors propose a Strecker reaction of the thiazolidine's open Schiff base giving **HSCH2CHO + CO2**
   (p. 84), i.e. mercaptoacetaldehyde (MAA). Cooked-meat aroma develops in the vials.
3. **The liberated thiol accumulates and persists** for 25 h at 55 °C — no sign of MAA collapsing to H2S on that
   timescale at that temperature.
4. For simple thiols (ME, GSH), the rate is set by a **spontaneous, pH-independent conversion of DH** to an activated
   intermediate that then takes up 2 thiols fast: k1 = 0.029 (ME) and 0.020 (GSH) h⁻¹ at 55 °C, and glycine
   accelerates it (k2 ≈ 1.3 M⁻¹) via a Schiff base (NaBH3CN suppresses it ~20×).

## 2. Methods as they matter to a model

- **DH.** Prepared after Madson & Feather 1981 (bis(benzoylhydrazone) route), standardised by conversion to
  metasaccharinic acid (Anet 1961).
- **Charge.** DH 5-50 mM, thiol 5-50 mM, glycine 0-1.0 M; "pH 5·5 acetate buffer" with [CH3COOH] + [CH3COO⁻] = 50 mM,
  **EDTA 0.5 mM** (essential: 10 µM MnCl2 greatly accelerated thiol loss — metal-catalysed autoxidation).
  NaBH3CN runs (50-75 mM) in **pH 6.9 phosphate**, [KH2PO4] = [Na2HPO4] = 25 mM.
- **Vessel / T.** ~2 mL in 2 mL screw-cap vials, water bath **55.0 ± 0.1 °C**; vials withdrawn and analysed.
  10 °C runs: reactants equilibrated, mixed, **sampled at 3 min intervals**, diluted at 10 °C and analysed at once.
- **Analysis.** Ellman (DTNB at pH 8, 412 nm), calibrated against iodine-standardised thiol (~0.05 M). Free thiol:
  acidify (1 mL 4 M HCl per 10 mL) and titrate with 10 mM iodine; total thiol: saturate with NaHCO3, titrate.
  **The Ellman assay reads free + reversibly bound thiol for ME/GSH** (hemithioacetals dissociate), but the
  **cysteine adduct "is not decomposed"** by it (p. 83). Replicates within ± 2 %.
- **TLC.** Silica gel G, butanol : acetic acid : water 4 : 1 : 1, DTNB spray.
- **Observable everywhere:** thiol concentration only. **DH, the adduct and the liberated thiol are never measured
  directly**; the liberated thiol's identity is not established.

## 3. Tables re-typed and figures digitised

### Table 1 (p. 76). k1 and k2 in k = k1(1 + k2[glycine]) for DH + ME / GSH at 55 °C (pH 5.5 acetate)

| thiol | k1 / h⁻¹ | k2 / M⁻¹ |
|---|---|---|
| mercaptoethanol | 0.029 | 1.31 |
| glutathione | 0.020 | 1.35 |

Supporting text (p. 75): first order in DH, zero order in thiol over 10-50 mM, at [gly] = 0 or 0.5 M. pH 3.85-6.00 at
[DH] = [ME] = 25 mM, [gly] = 0.5 M: initial rate **(1.06 ± 0.06) × 10⁻³ M h⁻¹** (n = 5) — **zero order in H⁺**.
At 0.5 M glycine: **k = 0.044 ± 0.003 h⁻¹ (ME, n = 15), 0.032 ± 0.002 h⁻¹ (GSH, n = 6)**; in D2O 0.045 and 0.033 h⁻¹
(no solvent isotope effect). Rate equation (p. 76): −d[SH]/dt = k′[DH](1 + k2[glycine]). NaBH3CN 50 mM cuts the
initial ME rate (DH = ME = 25 mM, gly 0.5 M) to ~1/20.

**Overall reaction (p. 76, Fig. 2):** DH 15.3 mM, ME 49.9 mM, gly 0.5 M, c∞ = 19.8 mM: first order over ≥ 90 % with
**k = 0.0106 h⁻¹**, about one quarter of the initial-rate constant. Stoichiometry DH : thiol **1 : 1.97 (ME), 1 : 2.03
(GSH)**. Derived (p. 81): first-order loss of DH **0.022 (ME) and 0.016 (GSH) h⁻¹** at 0.5 M glycine.

**For comparison (p. 72, recalculated from Wedzicha & Kaban 1986):** DH + S(IV), k = 5.32 × 10⁻³(1 + 4.45[glycine]) h⁻¹,
pH 5.5, 55 °C.

**Reversible binding (p. 74):** DH 5.2 + ME 19.5 mM + gly 0.5 M, 21 h: free thiol (iodine) 8.4 mM, total (iodometric,
raised pH) 16.9 mM, Ellman 16.6 mM. GSH, 24 h: 8.8 / 15.8 / 15.6 mM. Reversibly bound thiol 8.5 (ME) and 7.0 (GSH) mM
— **more than [DH]0**.

### Table 2 (p. 78). Cysteine bound to DH after 1 h; [Cys] = 24.9 mM, [gly] = 0.5 M, 50 mM acetate pH 5.5, 55 °C

| [DH]0 (mM) | 9.7 | 24.6 | 38.8 | 48.5 |
|---|---|---|---|---|
| [cysteine]bd (mM) | 8.5 | 14.9 | 18.7 | 20.8 |
| [cysteine]bd / [DH]0 | 0.88 | 0.61 | 0.48 | 0.43 |

Text (p. 77-78): extrapolating the 9.7 mM run to t = 0 gives 9.3 mM bound → **DH : cysteine 1 : 0.96**.

### Table 3 (p. 79). Initial rate of cysteine + DH vs [cysteine]; [DH] = 24.4 mM, 50 mM acetate pH 5.5, **10 °C**, no glycine

| [cysteine] (mM) | 10 | 15 | 25 | 40 | 50 |
|---|---|---|---|---|---|
| 10⁴ rate / M min⁻¹ | 3.90 | 3.94 | 3.79 | 4.08 | 3.46 |

### Fig. 3 [D] — total thiol (mM) vs time (h); [Cys]0 = 24.9 mM, gly 0.5 M, acetate pH 5.5, 55 °C

Read error ± 0.2 mM, ± 0.1 h. First samples are at ~0.1-0.2 h; the t = 0 value (24.9) is the charge, not plotted.

| [DH]0 | ~0.2 h | ~0.6 h | ~1.0 h | ~1.6 h | ~2.6 h | ~4.1 h | ~5.1 h | ~6.1 h | ~6.6 h | 25 h |
|---|---|---|---|---|---|---|---|---|---|---|
| 9.7 mM (△) | 19.1 | — | 16.4 | — | 16.5 | 16.8 | — | — | 17.7 | 15.8 |
| 24.6 mM (○) | 17.1 | 10.1 | 9.9 / 10.2 | 10.3 | 11.6 / 11.3 | 13.1 / 12.7 | 13.7 | 14.7 | 14.7 | 16.2 |
| 38.8 mM (▽) | 14.3 | — | 6.2 | — | 8.2 | 10.6 | — | — | 12.1 | 15.3 |
| 48.5 mM (□) | 12.9 | — | 4.1 | — | 7.1 | 9.4 | — | — | 11.1 | 11.8 |

(Duplicate markers at one time are printed as "a / b".) The 1 h minima reproduce Table 2 to ≤ 0.1 mM.

### Fig. 4 [D] — total thiol (mM) vs time (min); [Cys] = [DH] = 25 mM, 50 mM phosphate **pH 6.9**, 55 °C, no glycine

Read error ± 0.1 mM, ± 2 min.

| t (min) | 22 | 43 | 62 | 111 | 191 | 230 | 310 | 360 | 388 |
|---|---|---|---|---|---|---|---|---|---|
| no NaBH3CN (○) | 6.6 | 6.0 | 5.9 | 6.3 | 7.6 | 8.3 | 10.0 | 11.2 | 11.9 |
| NaBH3CN 75 mM (△) | 9.2 | 9.7 | 10.1 | 10.0 | 10.6 | 11.0 | 11.6 | 11.9 | 12.2 |

The authors' extrapolated lines: intercepts **4.4 mM (○) and 9.15 mM (△)** [D]; slopes **0.0169 and 0.0075 mM min⁻¹**
[D] (ratio 2.25; text: "2·2"); initial loss reduced "by some 25%" (20.6 → 15.85 mM bound [D], −23 %).

### Arithmetic (mine)

**1. Binding at 10 °C as a first-order constant.** Mean of Table 3 = 3.83 × 10⁻⁴ M min⁻¹ at 24.4 mM DH, zero order in
cysteine → if first order in DH (**not tested at 10 °C**), **k_bind ≈ 0.016 min⁻¹ (0.94 h⁻¹) at 10 °C, pH 5.5**. For
scale: the ME/GSH rate-determining DH conversion at **55 °C** is 0.020-0.029 h⁻¹ — cysteine's binding at 10 °C is
**~30-50× faster** than that, so cysteine does not wait for the "activated DH" that ME needs. Zero order in cysteine at
10 mM means the limiting step is a change of DH itself (plausibly ring-opening/dehydration of the cyclic/hydrated
forms to the free carbonyl) — the authors' "slow rearrangement of DH, e.g. involving keto-enol forms" (p. 84).

**2. Apparent binding quotient at 1 h, 55 °C** (Table 2; Q = bound / (free Cys × free DH), free DH = [DH]0 − bound):
DH 9.7 → **430 L/mol**; 24.6 → 154; 38.8 → 150; 48.5 → 183. **Not constant**; the 1 h point is past the onset of
thiol liberation, glycine (0.5 M) competes for DH, and the 9.7 mM arm is dominated by the small free-DH difference.
Read as **"K ≥ ~150 L/mol at 55 °C, pH 5.5"** at most — class **derived_from_table, low confidence**. Compare Cao
2017's TTCA (pentose + Cys) K_app ≈ 60 L/mol at 90 °C, pH ~4.5: the 3-deoxyosone binds cysteine at least as strongly
as the parent sugar.

**3. Liberation of the new thiol, 55 °C.** Taking bound pool B(t) = 24.9 − thiol(t) and the 1 h minimum as start:

| arm | liberated 1 → 6.6 h (mM) | mean rate (mM h⁻¹) | pseudo-first-order on B (h⁻¹) |
|---|---|---|---|
| DH 9.7 mM, pH 5.5 + gly | 1.3 | 0.23 | 0.030 |
| DH 24.6 mM, pH 5.5 + gly | 4.7 | 0.84 | 0.068 |
| DH 38.8 mM, pH 5.5 + gly | 5.9 | 1.05 | 0.068 |
| DH 48.5 mM, pH 5.5 + gly | 7.0 | 1.25 (2.0 over 1-2.5 h) | 0.073 |
| Cys = DH = 25 mM, pH 6.9, no gly (Fig. 4) | 6.0 (62 → 388 min) | 1.02 (extrapolated slope) | ~0.05 |

**~0.07 h⁻¹ (1.2 × 10⁻³ min⁻¹) at 55 °C** when DH ≥ cysteine; lower (0.03 h⁻¹) when DH is limiting, i.e. the
liberation step is faster with free DH in excess — consistent with a second DH-dependent step. **pH 5.5 → 6.9 changes the
rate by < 1.5×** (with the glycine confound). At 25 h three of four arms plateau near 15-16 mM, i.e. **~9 mM of the
cysteine is never liberated as thiol**, either held or converted to a non-thiol.

**4. NaBH3CN.** Reducing the Schiff base would, if the Strecker requires ring-opening to the Schiff base, slow liberation;
it halves it (ratio 2.2). This is the authors' evidence for a Schiff-base-dependent liberation step — weak (the reagent
also reduces free DH carbonyls at pH 6.9 slowly).

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | evidence class |
|---|---|---|---|---|---|
| DH → activated DH (thiol-independent RDS) | k1 | **0.029 (ME), 0.020 (GSH)** | h⁻¹ | 55 °C, pH 5.5 acetate, no gly | **measured** (Table 1) |
| glycine catalysis of the same | k2 | **1.31 (ME), 1.35 (GSH)** | M⁻¹ | 55 °C, pH 5.5 | measured (Table 1) |
| pH dependence of the RDS | order in H⁺ | **0** over pH 3.85-6.00 | — | 55 °C, gly 0.5 M | measured |
| DH + 2 RSH overall | first-order k | 0.0106 | h⁻¹ | DH 15.3, ME 49.9, gly 0.5 M | measured |
| stoichiometry | DH : RSH | 1 : 1.97 (ME), 1 : 2.03 (GSH); **1 : 0.96 (Cys, initial)** | — | 55 °C | measured |
| Cys + DH → adduct (thiazolidine) | initial rate | **3.8 × 10⁻⁴** | M min⁻¹ | 10 °C, pH 5.5, DH 24.4 mM, **zero order in Cys** 10-50 mM | measured (Table 3) |
| same | first-order k on DH | **~0.016** | min⁻¹ | 10 °C | derived (DH order assumed) |
| Cys–DH adduct extent | bound / DH0 at 1 h | 0.88 / 0.61 / 0.48 / 0.43 | — | 55 °C, gly 0.5 M, Cys 24.9 | measured (Table 2) |
| adduct → new thiol (proposed MAA) | pseudo-first-order on bound | **~0.07** (0.03-0.073) | h⁻¹ | 55 °C, pH 5.5 + gly; ~0.05 at pH 6.9 | **derived_from_digitised** |
| same, NaBH3CN effect | rate ratio | 2.2 | — | pH 6.9 | measured (text) |
| thiol vs t, cysteine arms | time series | Fig. 3, Fig. 4 tables | mM | 55 °C | digitised |

No activation energy for any step (55 °C for kinetics, 10 °C for binding only).

## 5. What this can and cannot anchor

**Can:**
1. **A missing edge: Cys + 3-deoxyosone ⇌ adduct**, fast and 1:1. In the engine, TDP cannot bind cysteine. Zero order in
   cysteine with a DH-intrinsic RDS → a first-order-in-osone formation step whose rate is set by the osone, ~0.016 min⁻¹ at
   10 °C, pH 5.5 (hexose analogue).
2. **A rate for the adduct's slow transformation** that liberates a non-cysteine thiol: ~0.07 h⁻¹ at 55 °C, weakly
   pH-dependent between 5.5 and 6.9. If that thiol is MAA (the authors' proposal; unconfirmed), it is the **first rate for
   a Strecker-type MAA release from a cysteine–deoxyosone adduct**, and the 25 h plateau says MAA (or whatever the thiol
   is) **persists as a thiol for a day at 55 °C**.
3. **pH-independence of the deoxyosone's RDS** (3.85-6.00) — relevant to how the engine weights deoxyosone reactions
   by pH.

**Cannot:**
- No pentose deoxyosone; no temperature series; no H2S measured; no product identified (MAA proposed, never isolated;
  the authors say they are "continuing with the synthesis" of MAA to test it).
- Glycine (0.5 M) is present in every 55 °C cysteine arm except Fig. 4; the bound/released numbers are glycine-confounded.
- Ellman readouts conflate cysteine, MAA, and any other thiol: the "liberated thiol" is a sum.

## 6. Flags

1. **The released thiol is not cysteine** (TLC Rf 0.12 vs 0.27). Any model reading "slow release of cysteine from the
   thiazolidine" from this paper misreads it; the paper says a different thiol forms.
2. **DH ≠ TDP.** Hexose 3-deoxyosone; the engine's TDP is the pentose analogue. Transfer by analogy only.
3. **Zero order in cysteine** means the formation rate does not scale with cysteine; a second-order mass-action step
   (k[Cys][osone]) fitted to this would be wrong in form at these concentrations.
4. **Apparent K at 1 h is not an equilibrium constant** (arithmetic 2).
5. **55 °C ≪ 100-145 °C** of the fit pots; no Ea. Extrapolation needs a declared Ea band.
6. The reversibly bound ME exceeds [DH]0 (p. 74): DH or its products carry more than one carbonyl/thiol site.
7. Registry: paper absent from `data/keys/papers.yml`; DH (3-deoxyglucosone), MAA, the adducts are not engine species.

## Use in the model

- **Gap (b), de Roos path.** Direct evidence that cysteine binds a **3-deoxyosone** 1:1, fast (minutes at 10 °C,
  pH 5.5), with a formation rate set by the osone (zero order in Cys), and that the adduct is not inert: it converts over
  hours at 55 °C (~0.07 h⁻¹) liberating a thiol. For B59 this argues for a `TDP + Cys ⇌ TDP-thiazolidine` edge alongside
  `PENT + Cys ⇌ TTCA`; at present TDP has no cysteine reaction at all. Use the binding as a **declared assumption** (hexose
  analogue, 10-55 °C) — not a fit row.
- **Gap (a), Strecker sulfur.** The best evidence in the corpus that a Strecker-type step **from the cysteine-osone
  adduct** produces a **thiol** (proposed MAA) rather than H2S directly, at 55 °C; the thiol persists 25 h. That supports
  routing the Strecker sulfur to an explicit MAA pool (not to H2S instantaneously, and not to inert `FRAG_S` unconditionally),
  with MAA's further fate (H2S + acetaldehyde, reductant-dependent per `kobayasi1965_extraction.md`) as a separate,
  temperature-dependent step. **No H2S was measured**, so the H2S split remains unmeasured.
- **Contradiction with current assumptions:** the engine's only cysteine-dicarbonyl Strecker is with the **1-deoxyosone**
  (`r_dpo_ddp`, second order). This paper's cysteine chemistry is with the **3-deoxyosone**, and its kinetics are
  binding-first, osone-limited, then slow first-order conversion — not a single bimolecular Strecker step.
- **Caveats:** hexose osone; 55 °C (kinetics) and 10 °C (binding) only; glycine 0.5 M in most cysteine arms; thiol is
  the only observable; product identities proposed, not shown.

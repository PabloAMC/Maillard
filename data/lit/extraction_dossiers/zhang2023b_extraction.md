# Zhang, Wang & Cao 2023 — EXTRACTION
### *"Factors influencing on the formation of dimethyl disulfide and dimethyl trisulfide in model systems"*
### Methionine 2.5 mg/mL (16.8 mM) in 0.01 M phosphate pH 4.9, sealed vials, **115 °C**, 30-150 min; ± lipids, ± 16 amino acids (1:1 w/w), ± reducing sugar (Xyl/Glc/Fru), ± cysteine or purified Cys-sugar Maillard intermediates (TTCA, ARP, HRP), initial pH 3-8. HS-SPME + GC-SIM with isopropyl disulfide internal standard and external standard curves (Table S2, not on disk). **Met alone: Met 14.9 -> 9.2 mM over 150 min; DMDS, DMTS, methional, MeSH all rise (sub-ng/mL to ~1 ng/mL). Cysteine or Cys-MRPs suppress DMDS/DMTS 2 to >100×; sugars + amino acids raise them up to ~300 ng/mL DMDS at pH 8.** Nearly every number is read off a figure; the absolute levels disagree between figures by up to ~700× (flag 1). Read 2026-10-28.

**Full read of `data/articles/zhang2023b.pdf`** (9 pp., typeset). Text layer clean; Table 1 transcribed from
`pdftotext -layout`. **Figs. 1a, 2a-j, 3a-h, 5a-d read off the embedded 300-dpi rasters (pp. 4-7)** by me:
Fig. 1a by pixel calibration against both axes (±0.005 ng/mL, ±0.3 mM); bar charts ±3 % of axis span.
Fig. 1b (correlation matrix) transcribed; Fig. 4 (inhibition rates, derived from Fig. 3) not digitised.
Supplementary Tables S1-S3 and Figs. S1-S2 (Met residuals, pH values, storage volatiles) not on disk.

**Provenance:** **[M]** printed number · **[F]** read off a figure by me · **[Q]** qualitative · **[Z]** derived by me · **[!]** flag.

---

## 0. IDENTITY — matches (with a repo-side mapping conflict)

| field | value |
|---|---|
| Authors | Zeyu Zhang, Bei Wang*, Yanping Cao* |
| Title | "Factors influencing on the formation of dimethyl disulfide and dimethyl trisulfide in model systems" (first page matches the requested paper) |
| Venue | *Food Research International* **172 (2023) 113200** (printed; Crossref agrees, issue Oct 2023) |
| DOI | **10.1016/j.foodres.2023.113200** (printed) |
| Dates | received 4 Jan 2023; revised 22 Jun 2023; accepted 27 Jun 2023; online 28 Jun 2023 |
| Affiliation | Beijing Technology and Business University (BTBU) |
| Repo note | [!] `docs/guides/CITATION_CONTEXT_2026-10.md` maps `zhang2023b.pdf` to this paper (correct), but `docs/guides/LITERATURE_SEARCH_2026-10.md` row 83 maps `zhang2023b.pdf` to a different paper (Zhang 2023, JAFC, lysine-xylose Amadori, 10.1021/acs.jafc.3c00388). The PDF on disk is this DMDS/DMTS paper. Not changed here |
| Companions | Pan 2025 (`pan2025_extraction.md`), Deng 2022 (`deng2022_extraction.md`), Yao 2025 (`yao2025_extraction.md`), Chin & Lindsay 1994 — the methionine-lane sources named in `src/kinetic_core/parameters_methionine.py` |

---

## 1. THE ONE-PARAGRAPH ANSWER

Methionine heated alone at 115 °C in pH 4.9 phosphate loses ~38 % in 150 min (14.9 -> 9.2 mM [F]; [Z]
apparent first-order k ≈ 3.2 × 10⁻³ min⁻¹) while methional, methanethiol, DMDS and DMTS rise together
(pairwise r 0.81-0.99, all r ≈ −0.9 with Met [M]). Adding **cysteine, or Maillard products/intermediates
made from cysteine**, accelerates Met loss yet **suppresses DMDS and DMTS** at pH 3-7 (and lowers pH);
adding a sugar with a non-sulfur amino acid raises DMDS 10-100× (e.g. Met-Xyl-Lys ~100, Met-Fru-Glu ~79
ng/mL DMDS [F]); raising the initial pH to 8 raises DMDS 5-50× in every sugar system (Glc-Met ~300 ng/mL
[F]). Cys-xylose TTCA leaves methanethiol high (~500-790 ng/mL) but DMDS low (≤ 0.5 ng/mL) [F], i.e. it
blocks the MeSH -> DMDS oxidation step. **But the absolute concentrations are not internally consistent:**
the same Met-alone 60-min pot reads DMDS 0.10 (Fig. 1a) vs 0.52-0.76 ng/mL (Figs. 2, 3, 5), methional 0.024
vs 17.7 ng/mL, and DMTS shows a ~1.1 ng/mL floor even in unheated Met (flags 1-2). Use the directions and
within-figure ratios, not the absolute ng/mL.

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Reagents | Met and 16 amino acids, Xyl, Glc, Fru (≥ 98 %); "cystine (Cys)" in the reagent list [!] (cysteine elsewhere); standards: methional, DMDS, DMTS (TCI), MeSH (Macklin), isopropyl disulfide (IPDS, internal standard) |
| MRIs | Cys 0.2 mol/L + sugar 1:2 molar, 50 mL water, pH 7.4 ± 0.1 (6 M NaOH), **115 °C 60 min**; Dowex 50WX4 (H+), water wash, 0.2 M NH4OH elution; lyophilised. Purities: XC TTCA (2-threityl-thiazolidine-4-carboxylic acid) ≥ 99 %, GC ARP ≥ 99 %, FC HRP ≥ 96 % (UPLC-MS, ¹³C NMR in Table S1) |
| Base pot | **Met 2.5 mg/mL** in sodium phosphate **pH 4.9, 0.01 M**; 30-mL sealed pressure vials; **115 °C, 60 min**; cooled under running water |
| Amino-acid series | Met : amino acid **1:1 (w/w)** (16 amino acids) |
| MRP series | "Met/Cys/various reducing sugars ... in a ratio of **1:1:7.5**" (pre-experiment optimum); Met : Cys-MRI 1:1 (w/w) in water; initial pH adjusted to **3, 4, 5, 6, 7, 8** with 0.01 M NaOH / H3PO4 |
| Time series | Met alone 115 °C for **30, 60, 90, 120, 150 min**; then mixed with Xyl/Cys ("Xyl:Cys:Met = 1:1:7.5, w/w") or XC TTCA (Met : TTCA 1:1 w/w) and heated **another 60 min at 115 °C** (Fig. 5) |
| Lipids | Met 2.5 mg/mL + **0.94 % (w/w)** tallow or soybean oil, same buffer and heating |
| Storage | Xyl-Cys-Met in several proportions (Table S3) at pH 4.9, **50 °C for 1, 28, 56 days** |
| Amino acids | AccQ-Tag pre-column derivatisation, HPLC-fluorescence |
| Volatiles | 10 mL solution + 1 µL IPDS ("0.0943 ng/mL") [!] in 40-mL vial; 45 °C 30 min equilibration + 30 min DVB/CAR/PDMS 50/35 µm; DB-WAX 60 m; full scan for ID (NIST 14); **SIM for MTP (m/z 48, 104), MeSH (47, 48), DMDS (79, 94), DMTS (79, 126), IPDS (108, 150)**; standard curves of area ratio vs IS (Table S2, not on disk). Triplicate |
| pH | hand-held meter after reaction |

---

## 3. TABLES AND FIGURE READINGS

### Fig. 1a. Methionine alone, pH 4.9, 115 °C [F] (Met on right axis, mM; odorants on left axis, ng/mL)

| time (min) | Met (mM) | methional (MTP) | MeSH | DMDS | DMTS |
|---|---|---|---|---|---|
| 0 | 14.9 | 0.009 | 0.034 | 0.050 | 0.008 |
| 30 | 13.6 | 0.015 | 0.051 | 0.069 | 0.015 |
| 60 | 12.4 | 0.024 | 0.104 | 0.103 | 0.034 |
| 90 | 11.0 | 0.025 | 0.105 | 0.136 | 0.044 |
| 120 | 10.1 | 0.024 | 0.134 | 0.194 | 0.056 |
| 150 | 9.2 | 0.030 | 0.140 | 0.327 | 0.090 |

[Z] Nominal Met 2.5 mg/mL = 16.8 mM (MW 149.21); t = 0 reads 14.9 mM (89 %). Apparent first-order loss
ln(14.9/9.2)/150 min = **3.2 × 10⁻³ min⁻¹** (5.4 × 10⁻⁵ s⁻¹) at 115 °C, pH 4.9; a 0-150 min straight line
fits the bars about as well (−0.038 mM/min). Error bars on Met are up to ±1.5 mM.

### Fig. 1b. Pearson correlations (n = time points; * p ≤ 0.05, ** p ≤ 0.01) [M]

| | MTP | MeSH | DMDS | DMTS | Met |
|---|---|---|---|---|---|
| MTP | – | 0.95 ** | 0.81 * | 0.89 * | −0.93 ** |
| MeSH | | – | 0.85 * | 0.92 * | −0.88 * |
| DMDS | | | – | 0.99 ** | −0.86 * |
| DMTS | | | | – | −0.92 * |

### Fig. 2a-b. Lipids (Met 2.5 mg/mL, 0.94 % w/w lipid, pH 4.9, 115 °C 60 min) [F], ng/mL

| group | DMDS | DMTS |
|---|---|---|
| Met-0h (unheated) | 0.14 (d) | 1.11 (c) |
| Met-1h | 0.76 (c) | 1.30 (b) |
| + tallow | 0.96 (b) | 1.29 (b) |
| + soybean oil | 1.44 (a) | 1.37 (a) |

### Fig. 2c-d. Free amino acids (Met : aa 1:1 w/w, 115 °C 60 min) [F], ng/mL (text: Met-Arg "reaching 3.32 ng/mL and 1.80 ng/mL" [M])

| group | DMDS | DMTS | | group | DMDS | DMTS |
|---|---|---|---|---|---|---|
| Met-0h | 0.14 | 1.11 | | Met-Ala | 0.58 | 1.28 |
| Met-1h | 0.76 | 1.30 | | Met-Pro | 0.85 | 1.23 |
| Met-Asp | 1.37 | 1.26 | | **Met-Cys** | **0.37** | **1.18** |
| Met-Ser | 0.64 | 1.27 | | Met-Tyr | 0.33 | 1.21 |
| Met-Glu | 0.28 | 1.15 | | Met-Val | 1.14 | 1.29 |
| Met-Gly | 0.72 | 1.24 | | Met-Lys | 0.50 | 1.19 |
| Met-His | 0.81 | 1.21 | | Met-Ile | 1.19 | 1.30 |
| **Met-Arg** | **3.32** | **1.80** | | Met-Leu | 0.78 | 1.36 |
| Met-Thr | 0.75 | 1.27 | | Met-Phe | 1.23 | 1.26 |

### Fig. 2e-j. Met + sugar + amino acid (ratio 1:1:7.5 as printed; 115 °C 60 min) [F], ng/mL

"Met-sugar-1h" = Met + sugar without a third amino acid.

| third component | Xyl: DMDS | Xyl: DMTS | Glc: DMDS | Glc: DMTS | Fru: DMDS | Fru: DMTS |
|---|---|---|---|---|---|---|
| Met-0h | ~0.1 | 1.1 | 0.15 | 1.11 | ~0.1 | 1.1 |
| Met-1h | ~0.8 | 1.3 | 0.5 | 1.15 | 0.5 | 1.2 |
| (none) Met-sugar-1h | 12 | 3.1 | 2.05 | 1.32 | 15 | 4.5 |
| Asp | 21 | 6.4 | 1.95 | 1.35 | 21 | 6.7 |
| Ser | 24 | 3.8 | 1.5 | 1.27 | 9 | 2.8 |
| Glu | 50 | 15.4 | 6.6 | 2.0 | 79 | 29.6 |
| Gly | 49 | 11.0 | 1.95 | 1.36 | 18 | 6.1 |
| His | 22 | 27.7 | 1.65 | 1.45 | 6.5 | 5.0 |
| Arg | 84 | 29.0 | 4.85 | 2.0 | 36 | 18.4 |
| Thr | 38 | 8.4 | 2.4 | 1.43 | 17 | 6.7 |
| Ala | 25 | 5.8 | 5.1 | 1.44 | 19.5 | 7.0 |
| Pro | 24 | 5.5 | 6.55 | 1.67 | 8.5 | 3.1 |
| **Cys** | **1.5** | **1.2** | **0.3** | **1.13** | **0.5** | **1.1** |
| Tyr | 18 | 7.9 | 1.75 | 1.74 | 4.5 | 6.1 |
| Val | 35 | 9.2 | 7.1 | 1.55 | 13 | 4.2 |
| Lys | 100 | 29.3 | 13.8 | 2.93 | 40 | 18.5 |
| Ile | 25 | 5.2 | 2.15 | 1.36 | 11 | 3.1 |
| Leu | 61 | 26.7 | 6.7 | 5.4 | 43 | 26.3 |
| Phe | 32 | 7.9 | 6.0 | 1.57 | 28 | 10.3 |

### Fig. 3. Effect of initial pH (115 °C 60 min) [F], ng/mL

| system | pH 3 | pH 4 | pH 5 | pH 6 | pH 7 | pH 8 |
|---|---|---|---|---|---|---|
| **DMDS** Met | 0.36 | 0.29 | 0.52 | 0.60 | 0.67 | 3.05 |
| DMDS Met-Cys | 0.14 | 0.15 | 0.14 | 0.17 | 0.13 | 0.55 |
| DMDS Xyl-Met | 4 | 5 | 11 | 18 | 33 | 140 |
| DMDS Xyl-Cys-Met | ~0.5 | ~0.5 | ~0.5 | ~0.5 | ~0.5 | 9 |
| DMDS XC TTCA-Met | ~0.5 | ~1 | ~1.5 | ~2 | ~2.5 | 3.5 |
| DMDS Glc-Met | 4 | 3 | 7 | 6 | 28 | 302 |
| DMDS Glc-Cys-Met, GC ARP-Met | < 1 | < 1 | < 1 | < 1 | < 1 | < 1 |
| DMDS Fru-Met | 5 | 5 | 5.5 | 14 | 18 | 176 |
| DMDS Fru-Cys-Met, FC HRP-Met | < 1 | < 1 | < 1 | < 1 | < 1 | < 1 |
| **DMTS** Met | 1.24 | 1.19 | 1.23 | 1.20 | 1.18 | 1.32 |
| DMTS Met-Cys | 1.14 | 1.15 | 1.14 | 1.16 | 1.15 | 1.97 |
| DMTS Xyl-Met | 1.65 | 1.88 | 2.58 | 3.10 | 4.75 | 6.80 |
| DMTS Xyl-Cys-Met | 1.18 | 1.17 | 1.15 | 1.16 | 1.22 | 1.70 |
| DMTS XC TTCA-Met | 1.17 | 1.24 | 1.40 | 1.62 | 1.85 | 4.55 |
| DMTS Glc-Met | 1.24 | 1.33 | 1.74 | 1.65 | 2.19 | 3.47 |
| DMTS Glc-Cys-Met | 1.12 | 1.12 | 1.13 | 1.12 | 1.15 | 1.14 |
| DMTS GC ARP-Met | 1.10 | 1.12 | 1.17 | 1.15 | 1.26 | 1.80 |
| DMTS Fru-Met | 1.58 | 1.85 | 1.73 | 3.40 | 3.26 | 2.87 |
| DMTS Fru-Cys-Met | 1.13 | 1.13 | 1.13 | 1.13 | 1.14 | 1.14 |
| DMTS FC HRP-Met | 1.13 | 1.14 | 1.16 | 1.17 | 1.22 | 1.53 |

### Table 1. "Regression analysis of the changes in the inhibition rate of DMDS and DMTS in the range of pH 3 to 7" [M] (y = inhibition rate, %; x = pH, as I read it)

| group | DMDS | R² | DMTS | R² |
|---|---|---|---|---|
| A: Cys-Met | y = 6.35x + 35.43 | 0.67 | y = -1.05x + 10.04 | 0.43 |
| B: Xyl-Cys-Met | y = 10.91x + 28.36 | 0.80 | y = 11.94x − 8.62 | 0.99 |
| C: Glc-Cys-Met | y = 2.64x + 79.99 | 0.55 | y = 9.35x − 19.13 | 0.91 |
| D: Fru-Cys-Met | y = 4.62x + 70.15 | 0.58 | y = 10.06x − 3.33 | 0.80 |

### Fig. 5. Met pre-heated 30-150 min at 115 °C, then + Xyl/Cys or + XC TTCA, then 60 min more at 115 °C [F], ng/mL ("Met" = control, no addition)

| Met pre-heat (min) | MTP Met | MTP Xyl-Cys-Met | MTP TTCA-Met | MeSH Met | MeSH Xyl-Cys-Met | MeSH TTCA-Met | DMDS Met | DMDS Xyl-Cys-Met | DMDS TTCA-Met | DMTS Met | DMTS Xyl-Cys-Met | DMTS TTCA-Met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 14.2 | 8.1 | 10.0 | 14 | 115 | 490 | 0.27 | 0.38 | 0.29 | 1.165 | 1.135 | 1.20 |
| 60 | 17.7 | 7.9 | 10.4 | 25 | 148 | 510 | 0.71 | 0.38 | 0.35 | 1.28 | 1.13 | 1.20 |
| 90 | 14.5 | 10.2 | 6.9 | 32 | 207 | 787 | 0.76 | 0.27 | 0.47 | 1.31 | 1.145 | 1.19 |
| 120 | 13.2 | 6.1 | 1.8 | 52 | 217 | 735 | 1.12 | 0.19 | 0.41 | 1.32 | 1.125 | 1.17 |
| 150 | 14.3 | 2.7 | 2.4 | 55 | 162 | 745 | 0.87 | 0.19 | 0.21 | 1.34 | 1.12 | 1.17 |

### Storage (Xyl-Cys-Met, pH 4.9, 50 °C, 56 days) [Q]
12 volatile sulfur compounds over 56 days (Table S3, not on disk). Day 1: Xyl+Cys gives only 2-furfurylthiol
and difurfuryl disulfide; Xyl+Met gives FFT [sic, "furfuryl mercaptan"], DMDS, DMTS and methional. FFT decays
to undetectable with time in all systems; furfuryl methyl sulfide and 2-thiophenecarboxaldehyde appear after
~1 month in Cys-rich mixtures; MTP, DMDS, DMTS rise with time in the ternary Met/Xyl/Cys mixtures but "far
less" than in Xyl/Met.

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Model context.** The engine's methionine lane (B22, `src/kinetic_core/parameters_methionine.py`) is fitted
on Pan 2025 (Met 0.27 mM in a 240 mM hexose pool, citrate pH 6.2, 100-140 °C); it carries methionine ->
methional (Strecker), methional -> MeSH (retro-Michael), 2 MeSH -> DMDS (apparent), and **refuses DMTS**
because sugar + methionine pots do not carry H2S.

**Can anchor (directional / ratio level):**
- **Cysteine suppresses DMDS and DMTS from methionine at pH 3-7** in every sugar system (Fig. 2c-j "Cys"
  column, Fig. 3 "-Cys-" rows: 2 to >100× lower DMDS), while accelerating Met loss. A sign constraint on any
  future coupling of the sulfur (cysteine) lane to the methionine lane.
- **TTCA (Cys-xylose intermediate) holds MeSH high but DMDS low** (Fig. 5: MeSH 490-790 vs DMDS ≤ 0.5
  ng/mL): consistent with a MeSH -> DMDS step that is oxidant-limited (radical-scavenging by Cys-MRPs) — the
  engine's DMDS constant is explicitly "apparent; no oxidant tracked".
- **pH dependence of DMDS:** monotone rise from pH 3-4 to 7 (Xyl-Met 4 -> 33 ng/mL) and a jump at pH 8
  (×4-11 over pH 7). Within-figure ratios are usable as a soft pH-shape check for the methionine lane.
- **A third (non-sulfur) amino acid changes DMDS from Met + sugar by 0.3-8×**: with xylose every one raises
  it (1.8× Ser/His to 8× Lys: 12 -> 100 ng/mL); with glucose/fructose some lower it (Fru-His 0.4×, Fru-Tyr
  0.3×, Glc-Ser 0.7×). Lys, Leu, Arg, Glu are consistently the strongest promoters. Direction only.
- **DMTS from methionine without cysteine** rises with heating in Met + sugar pots (Xyl-Met 1.65 -> 6.8
  ng/mL over pH 3-8; Fig. 2f up to 29 ng/mL). This is evidence against the strict form of
  `METHIONINE_NO_DMTS_REASON` ("needs hydrogen sulfide ... which the sugar path does not carry"), unless
  the DMTS signal is the ~1.1 ng/mL artefact (flag 2); above ~3 ng/mL it clearly is not.

**Cannot anchor:**
- **No absolute rate constants for the volatiles** (flag 1), no time course except Fig. 1a (suspect) and
  Fig. 5 (a two-stage protocol).
- **Met loss rate** (3.2 × 10⁻³ min⁻¹, Met alone) has no product balance: Met lost (~5.6 mM by 150 min)
  vs volatile S products (≤ ~10⁻³ mM even on the highest-reading figure, Fig. 5 MeSH; ~10⁻⁵ mM on Fig. 1a)
  differ by > 10³; the sink (sulfoxide, Met-derived non-volatiles,
  assay drift) is unmeasured. Not a fit target without more.
- Final pH values (Fig. S2) and Met residuals with additives (Fig. S1) are in the supplement, not on disk.

**Caveats / flags:**
1. [!] **Absolute levels inconsistent across figures** for nominally the same pot (Met alone, pH 4.9, 115
   °C, 60 min): DMDS 0.10 (Fig. 1a) vs 0.76 (Fig. 2a), 0.52 (Fig. 3a, pH 5) and 0.71 (Fig. 5c); DMTS 0.034
   vs 1.2-1.3; methional 0.024 (Fig. 1a) vs 17.7 (Fig. 5a); MeSH 0.10 vs 25. Fig. 1a is 7-700× below the
   others. Either an axis-unit error in Fig. 1a or different calibrations; I cannot resolve it from the
   paper. Figs. 2, 3, 5 agree with each other for DMDS/DMTS.
2. [!] **DMTS floor ≈ 1.1 ng/mL** in unheated Met (Met-0h 1.11) and in every suppressed system (1.10-1.18):
   a baseline (calibration intercept, blank or carry-over), not formation. Net heat-formed DMTS = printed −
   ~1.1. Net Met-alone DMTS at 60 min is ~0.2 ng/mL, and the Met-Arg "1.80" is ~0.7 net.
3. [!] **Internal-standard amount implausible:** 1 µL of IPDS at "0.0943 ng/mL" puts ~0.1 pg in the vial,
   far below SPME-GC-SIM detection; probably mg/mL or µg/µL. Absolute ng/mL depend on this.
4. [!] **Ratio ambiguity:** "Met/Cys/various reducing sugars ... 1:1:7.5" vs "Xyl:Cys:Met = 1:1:7.5,
   w/w" (opposite order) — whether sugar or Met is the 7.5 part is unclear; and the reagent list says
   "cystine (Cys)".
5. [!] Initial-pH series used 0.01 M NaOH/H3PO4 in 0.01 M phosphate (or water for the MRI mixtures): weak
   buffering; the final pH (Fig. S2) differs from the nominal and is lower with Cys (stated [Q]).
6. [!] HS-SPME on 10 mL aqueous at 45 °C with an IS calibration curve (Table S2, not available): no
   matrix-matched or isotope correction; MeSH by SPME is particularly lossy.
7. Figure read-offs are mine; bar values ±3 % of axis span, Fig. 1a ±0.005 ng/mL / ±0.3 mM.

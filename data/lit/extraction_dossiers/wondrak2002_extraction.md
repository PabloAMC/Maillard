# Wondrak 2002 — EXTRACTION (α-oxoaldehyde scavenging by thiols and by aminoguanidine; four **second-order** rate constants at 37 °C, pH 7.4, 10 mM phosphate, by HPLC disappearance of phenylglyoxal; plus a 96-well AGE-fluorescence screen of fifteen nucleophiles on a histone H1 / ADP-ribose glycation system)

### THE HEADLINE: this paper prints the corpus's first **measured second-order rate constant for a thiol reacting with an α-dicarbonyl** — L-cysteine + phenylglyoxal, **k₂ₙd = 0.63 ± 0.04 M⁻¹ s⁻¹ at 37 °C, pH 7.4** — and that constant belongs to the **ADDUCT branch, not the redox branch**. The product is a 2-acylthiazolidine (structure elucidated by ¹H/¹³C NMR and MALDI-TOF-MS for the D-penicillamine analogues). **No disulfide is formed, sought, or mentioned anywhere in this paper.** It therefore supplies no rate, no order and no barrier for the step `ch_redox_dicarbonyl` that `kinetic_core_b27_prereg.md` proposes.

**Source on disk:** `data/articles/wondrak2002.pdf` (13 pp., Biochemical Pharmacology 63 (2002) 361–373).
Read from the `pdftotext -layout` text layer, with **pages 4, 5, 6 and 9 re-read as page images** through the Read tool to recover the micro sign, which the text layer renders as `m` throughout (so "50 mM phenylglyoxal" in the text layer is **50 µM** in the printed page, and "600 mM methylglyoxal" is **600 µM**). Every concentration below is taken from the page image, not the text layer. Table 1 came through the image clean and is re-typed in full. **All four rate constants in this paper are printed in the caption of Fig. 6 (p. 369) and in the running text of §3.4 (p. 368) — there is no rate-constant table.** Figures 3, 4, 6, 7 and 8 are plots; only the caption numbers are typed here.
**Repo status before this dossier:** Wondrak 2002 is cited nowhere in `src/kinetic_core/parameters_sulfur.py`, nowhere in `data/lit/reaction_rules.yml`, and has no extraction dossier. The only mentions of "penicillamine" anywhere in the repository are two raw deep-research dumps under `data/research_corpus/`.

## 0. Identity

| field | value |
|---|---|
| Title | "Identification of α-dicarbonyl scavengers for cellular protection against carbonyl stress" |
| Authors | Georg T. Wondrak, Daniel Cervantes-Laurean, Michael J. Roberts, Jaber G. Qasem, Moonsun Kim, Elaine L. Jacobson, Myron K. Jacobson (corresponding, mjacobson@pharmacy.arizona.edu) — Department of Pharmacology and Toxicology, College of Pharmacy, Arizona Cancer Center, University of Arizona, 1515 North Campbell Avenue, Tucson, AZ 85724, USA |
| Venue | Biochemical Pharmacology **63** (2002) 361–373. Received 9 March 2001; accepted 26 July 2001 |
| DOI | **No DOI is printed in the PDF.** The front matter carries only `PII: S 0 0 0 6 - 2 9 5 2 ( 0 1 ) 0 0 9 1 5 - 7` (p. 361) and the ISSN line `0006-2952/02/$`. Recorded as absent per house rule; do not synthesise one. |
| Funding / interest | NIH CA43894, NS38496, and Niadyne Inc. **M.K.J. and E.L.J. are principals in Niadyne Inc.** (Acknowledgments, p. 371) |
| Computational content | **None.** There is no DFT, no molecular mechanics, no docking and no quantum chemistry anywhere in this paper, so the standing no-DFT policy excludes nothing here. |
| The α-dicarbonyls | **phenylglyoxal** (the kinetic probe, chosen because it is UV-active at 254 nm), **methylglyoxal** (the preparative and cell-culture α-dicarbonyl), **glyoxal** (cell culture only), **ADP-ribose** (the glycation driver, a phosphate-substituted pentose whose deoxypentosones are the in-situ α-dicarbonyl source) |
| The thiols | L-cysteine, L-cysteine-O-methylester, Nᵅ-acetyl-L-cysteine (NAC), glutathione (GSH), L-Cys-Gly, D,L-homocysteine, cysteamine, D-penicillamine and D,L-penicillamine (3,3-dimethyl-D-cysteine), L-ergothioneine, 2-thiobarbituric acid, thiourea |
| Naming | "α-oxoaldehyde" = α-dicarbonyl; "RCS" = reactive carbonyl species; "k₂ₙd" = the second-order constant in the paper's own units `/M s`; "AGE fluorescence" = λex 370 / λem 440 nm (plate reader 355/405 nm) |

## 1. Why it matters

**What the pre-registration needs, and what this paper does and does not give it.**
`results/validation/kinetic_core_b27_prereg.md` §3(b) proposes one new step, `ch_redox_dicarbonyl`,
**second order, consuming an α-dicarbonyl and delivering oxidising equivalents to `OX`**, on a fitted
constant and a declared barrier, taking from Whitfield & Mottram 1999 the mechanism that the pot's
own α-dicarbonyls are reduced to hydroxyalkanones while thiols go to disulfides. The pre-registration
states in its own §1 that Whitfield "prints no rate, no order and no barrier" for that step.

Wondrak 2002 is the nearest thing in this batch to an anchor for it, and it **anchors the wrong
branch**:

- **What it confirms.** The encounter is real, direct, fast, and **second order** — `−dc/dt =
  k₂ₙd[phenylglyoxal][scavenger]` (§2.8, p. 364), verified by the apparent first-order constant
  scaling with scavenger concentration and by the two reactant ratios (5:1 and 10:1) agreeing
  (§2.8, p. 365). The **order** the pre-registration assumes for `ch_redox_dicarbonyl` therefore has
  an independent precedent for this reactant pair. That is the whole of the support.
- **What it does not confirm.** Every one of the four constants is the rate of **adduct formation**.
  The paper isolates the products and elucidates them: methylglyoxal + D-penicillamine gives
  **2-acetyl-5,5-dimethyl-thiazolidine-4-carboxylic acid** ([M+Na]⁺ = 226 Da by MALDI-TOF-MS,
  ¹H NMR fully assigned, ¹³C carbonyls at δ 199 and δ 209 ppm, **no aldehyde proton**, §3.3, p. 368);
  phenylglyoxal + D-penicillamine gives **2-benzoyl-5,5-dimethyl-thiazolidine-4-carboxylic acid**
  (§2.7, p. 364 and §3.3, p. 368). Ring closure of the thiol sulfur onto the Schiff-base carbon
  (Fig. 5, p. 368). **The α-dicarbonyl is consumed, and the thiol is consumed, and the sulfur ends
  up inside a five-membered ring bonded to carbon.** No sulfur–sulfur bond is made. The words
  "disulfide", "cystine", "GSSG" and "oxidised thiol" do not occur in this paper at all, and no
  oxidised-thiol species is measured, assayed or discussed.
- **The direct consequence for B27.** The redox step the pre-registration wants to fit does not
  merely lack a constant — it has a **measured competitor for the same two reactants**. In a pot fed
  cysteine (`fed_nf_cys_MFT`, `whitfield_nf_cys_MFT`, the two pots the pre-registration §2 names as
  the "not defensible" zeros), any α-dicarbonyl the proposed `r_nf_dicarbonyl` source manufactures
  will be attacked by cysteine at a *measured* 0.63 M⁻¹ s⁻¹ (37 °C, pH 7.4) along a route that
  produces **no oxidising equivalent at all**. A fitted `ch_redox_dicarbonyl` constant, having no
  anchor, will be free to be set arbitrarily far above that; nothing in the objective would notice.
  **This paper is the reason a fitted redox constant needs a declared ceiling**, and it supplies a
  candidate ceiling in the right units at the wrong temperature (Flags 1).
- **The pre-registration's second structural fact is untouched.** "No reaction anywhere makes an
  α-dicarbonyl from the fed substrate" — Wondrak's α-dicarbonyls are either bought (phenylglyoxal,
  methylglyoxal, glyoxal) or generated from **ADP-ribose**, a pentose, by the deoxypentosone route
  already in the repository. Nothing here makes an α-dicarbonyl from norfuraneol, and nothing here
  bears on `r_nf_dicarbonyl`.

**Where the numbers land in the repository.** `parameters_sulfur.py` already carries `k_thioether`
(5.01e-4, order 2) for `R-SH + matrix electrophile site -> matrix-bound thioether`, whose own
dossier note records that Hofmann 2001 "prints no rate constant, no order and no barrier" and that
the 9.8e-4 /s is a derived first-order reading. **Wondrak's constants are the opposite case: they
are printed by the authors, in stated units, with a stated order and a stated standard deviation,
from a reaction the authors wrote down explicitly.** They are the cleanest measured thiol +
electrophile second-order constants in the corpus. They are also not for a matrix electrophile and
not at 140 °C (Flags 1, 2).

**What this paper does NOT give the repository:** any disulfide measurement; any temperature other
than 37 °C, hence **no activation energy for anything**; any pH other than 7.4; any constant for
methylglyoxal or glyoxal itself (the kinetics are phenylglyoxal only — Flags 3); any constant for
GSH, NAC, cysteamine or ergothioneine (screened but not kinetically characterised); any oxygen
dependence (the whole assay was designed to be oxygen-independent); any headspace, volatile or
aroma measurement; any food matrix.

## 2. Methods as they matter to a model

**The kinetic pot (§2.8, p. 364–365; §2.10, p. 365) — this is the pot the four constants come from.**

- **Reactants.** Phenylglyoxal at **50 µM**; carbonyl scavenger at **250 µM and 500 µM** (two runs
  per compound, i.e. 1:5 and 1:10 phenylglyoxal : scavenger). Scavengers kinetically characterised:
  **D-penicillamine, aminoguanidine, L-cysteine**. A fourth run substitutes **phenylacetaldehyde**
  (a mono-oxoaldehyde) for phenylglyoxal against D-penicillamine.
- **Buffer, pH, temperature.** **10 mM phosphate buffer, pH 7.4, 37 °C.** No ionic-strength
  statement, no chelator in the kinetic runs, no stated ionic strength beyond the 10 mM phosphate.
- **Atmosphere.** **Not stated for the kinetic runs.** (The *screening* assay, a different pot, was
  run both under argon with 5 mM DTPA and under air, and gave the same answer — Table 1.) There is
  no argon, no degassing and no oxygen exclusion described in §2.8.
- **Time.** Sampling every **20 s** for D-penicillamine; the plotted windows in Fig. 6 are **0–60 min**
  (panel A) and **0–12.5 min** (panel B).
- **How the reaction was quantified.** **Disappearance of phenylglyoxal by HPLC** — C4-Rainin-Microsorb
  4.6 mm × 250 mm, 300 Å, 5 µm; isocratic 20 % acetonitrile / 80 % water (0.1 % TFA); 1 mL/min;
  UV detection at **254 nm**; the measured quantity is the **area under the curve (AUC)** of the
  phenylglyoxal peak. Aliquots kept on dry ice until analysis.
- **How the constant was obtained.** Pseudo-first-order: `log(AUC) vs. time` has slope `k₁ₛₜ/2.303`;
  then **`k₂ₙd = k₁ₛₜ / [α-dicarbonyl scavenger]`** (§2.8, p. 365, printed exactly so). Values from
  the 5:1 and 10:1 ratios "were in good agreement" — **the paper does not print the two separately**,
  only the pooled result (Flags 4). Means ± SD of **four** measurements (Fig. 6 caption).
- **What is NOT measured in this pot.** Product formation. The kinetics follow only the loss of the
  α-dicarbonyl. **The thiol is never assayed, and no oxidised thiol species is ever looked for.**
  The adduct identity comes from a *separate, preparative* experiment at very different
  concentrations (below), not from the kinetic runs.

**The preparative adduct pots (§2.6 and §2.7, p. 364) — where the branch assignment is established.**

- *Methylglyoxal route:* D-penicillamine **350 mg (2.3 mmol)** in **50 mL of 0.20 M phosphate buffer,
  pH 7.4**, plus methylglyoxal (40 % in H₂O, 620 µL, **3.45 mmol**); stirred **37 °C for 24 h**;
  desalted on Amberchrome CG 71 ms, then QAE Sephadex 25 anion exchange on a water → 0.2 M NH₄HCO₃
  gradient; product characterised by ¹H NMR (D₂O, Varian Gemini-200) and MALDI-TOF-MS (Kratos Kompact
  Seq, positive ion, linear, α-cyano-4-hydroxycinnamic acid matrix).
- *Phenylglyoxal route:* phenylglyoxal **10 mM** final, D-penicillamine **20 mM** final, **50 mM
  KH₂PO₄ pH 7.4, room temperature**; **>90 % conversion of the phenylglyoxal peak into a single
  product peak in 40 min** by HPLC at 254 nm; product taken by preparative HPLC (gradient 0.1 % TFA →
  50 % acetonitrile) and analysed by ¹H NMR.
- **Reversibility, as the paper states it.** "The α-dicarbonyl adduct is **stable in water**, whereas
  the monocarbonyl-derived thiazolidine adducts form **reversibly** with subsequent release of the
  aldehyde" (Discussion, p. 371, citing ref. [44]). The paper's own §3.3 heading is "**Irreversible**
  trapping of α-dicarbonyl compounds by D-penicillamine". **No equilibrium constant and no reverse
  rate constant is measured** — the irreversibility claim rests on the cited literature and on the
  observed stability of the isolated product (Flags 5).

**The screening pot (§2.4, p. 363; Table 1, p. 366) — a different pot, and the source of the GSH and
cysteine level data.**

- **1.5 mg/mL histone H1** (isolated in-house from calf thymus, §2.3) + **1 mM ADP-ribose** in
  **50 mM KH₂PO₄, pH 7.4, 37 °C**, with **0.015 % NaN₃**; total volume **300 µL** per well on a
  96-well plate sealed with a watertight plastic sheet; test compound **1–10 mM**; **5 days**.
- **Atmosphere:** run **both** under argon with **5 mM DTPA** and under air, and the two agree
  (Table 1: 23 (1.1) vs 21 (0.9) at day 5). This is the paper's demonstration that its glycation
  system is **oxygen- and metal-independent**, which is the entire methodological point.
- **Readout:** AGE fluorescence, λex 355 / λem 405 nm on a Fluoroskan II (bandwidth 35 nm), against
  the 1 mL reference geometry read at λex 370 / λem 440 nm on a Hitachi F-2000. Two orthogonal
  controls: the **AGE–BSA test** for fluorescence quenchers (false positives), and **12 % SDS-PAGE
  with silver staining** for protein cross-linking.
- **AGE–BSA** was made separately: 1.6 g BSA + 3.0 g D-glucose in 10 mL of 0.5 M sodium phosphate
  pH 7.4 with 0.05 % NaN₃, filter-sterilised, **90 days at 37 °C in the dark** (§2.2).
- **What Table 1 measures is a fluorescence level, not a rate and not a thiol concentration.**
  Nothing in Table 1 is a kinetic quantity.

**The cell pot (§2.9, p. 365; Figs. 7–8, p. 369).** HaCat keratinocytes (2×10⁴/well) and CF3
fibroblasts (4×10⁴/well) on 6-well dishes, DMEM + 10 % FBS, 5 % CO₂, 37 °C; scavenger added **15 min
before** the α-dicarbonyl; **72 h** exposure; counted on a Coulter counter. Methylglyoxal and glyoxal
at 0 / 300 / 600 µM; scavengers at 1 mM; D-penicillamine dose–response at 0 / 0.5 / 1.0 mM. Means ±
SD of three samples. **Cell-count endpoints only; no intracellular thiol, no GSH, no GSSG assay.**

## 3. Tables re-typed

**This paper contains exactly one table.** Every rate constant is in a figure caption or in running
text and is re-typed under "Numbers printed in figure captions and running text" below.

### Table 1 (p. 366). "Screening of inhibitors of non-oxidative advanced glycation: AGE fluorescence on 96-well microtiter plate"

| Sample | AGE fluorescence ᵃ (day 0) | AGE fluorescence ᵃ (day 5) | AGE–BSA test ᵃ,ᵇ |
|---|---|---|---|
| Histone H1 blank | 1.1 (0.0) `[M]` | 1.5 (0.2) `[M]` | |
| **Complete reaction** | | | |
| – Under argon (+5 mM DTPA) | 1.1 (0.1) `[M]` | 23 (1.1) `[M]` | 32 `[M]` |
| – Under air | 1.1 (0.1) `[M]` | 21 (0.9) `[M]` | 30 `[M]` |
| **+Compound** | | | |
| *Aminoguanidine (mM)* 1 | 1.2 (0.1) `[M]` | 4.3 (0.2) `[M]` | |
| 5 | 1.2 (0.1) `[M]` | 2.0 (0.0) `[M]` | |
| 10 | 1.3 (0.1) `[M]` | 1.8 (0.0) `[M]` | 10 `[M]` |
| *Rutin (**µM**)* 200 | 1.0 (0.0) `[M]` | 3.0 (0.1) `[M]` | 5.7 `[M]` |
| *NADH (mM)* 5 | 52 (0.3) `[M]` | 31 (0.0) `[M]` | |
| *L-Cys-Gly (mM)* 1 | 1.3 (0.3) `[M]` | 17 (0.0) `[M]` | |
| 5 | 1.2 (0.2) `[M]` | 26 (0.4) `[M]` | |
| 10 | 1.2 (0.1) `[M]` | 40 (0.0) `[M]` | |
| *GSH (mM)* 1 | 1.3 (0.0) `[M]` | 19 (0.3) `[M]` | |
| 5 | 1.1 (0.3) `[M]` | 17 (0.3) `[M]` | |
| 10 | 1.1 (0.0) `[M]` | 14 (0.5) `[M]` | |
| *L-Cys (mM)* 1 | 1.2 (0.1) `[M]` | 14 (0.3) `[M]` | |
| 5 | 1.1 (0.1) `[M]` | 11 (0.2) `[M]` | |
| 10 | 1.2 (0.1) `[M]` | 9.0 (0.2) `[M]` | |
| *L-Cys-OMe (mM)* 1 | 1.2 (0.1) `[M]` | 15 (0.1) `[M]` | |
| 5 | 1.3 (0.1) `[M]` | 12 (0.8) `[M]` | |
| 10 | 1.2 (0.1) `[M]` | 9.0 (0.2) `[M]` | |
| *NAC (mM)* 1 | 1.2 (0.2) `[M]` | 11 (0.6) `[M]` | |
| 5 | 1.0 (0.1) `[M]` | 10 (0.4) `[M]` | |
| 10 | 1.1 (0.1) `[M]` | 9.4 (0.6) `[M]` | |
| *D,L-Homocysteine (mM)* 1 | 1.3 (0.0) `[M]` | 16 (1.1) `[M]` | |
| 5 | 1.2 (0.2) `[M]` | 22 (0.1) `[M]` | |
| 10 | 1.1 (0.1) `[M]` | 27 (1.7) `[M]` | |
| *Cysteamine (mM)* 1 | 1.0 (0.1) `[M]` | 10 (0.6) `[M]` | |
| 5 | 1.1 (0.1) `[M]` | 9.0 (0.7) `[M]` | |
| 10 | 1.2 (0.2) `[M]` | 6.0 (0.3) `[M]` | |
| *D,L-Penicillamine (mM)* 1 | 1.1 (0.2) `[M]` | 13 (0.0) `[M]` | |
| 5 | 1.1 (0.1) `[M]` | 2.7 (0.0) `[M]` | |
| 10 | 1.1 (0.1) `[M]` | 1.7 (0.0) `[M]` | 15 `[M]` |
| *D-Penicillamine (mM)* 1 | 1.2 (0.0) `[M]` | 11 (0.6) `[M]` | |
| 5 | 1.1 (0.1) `[M]` | 2.9 (0.3) `[M]` | |
| 10 | 1.1 (0.0) `[M]` | 1.4 (0.1) `[M]` | 12 `[M]` |
| *2-Thiobarbituric acid (mM)* 1 | 1.1 (0.1) `[M]` | 12 (0.4) `[M]` | |
| 5 | 1.2 (0.1) `[M]` | 5.5 (0.4) `[M]` | |
| 10 | 1.3 (0.1) `[M]` | 2.3 (0.0) `[M]` | 12 `[M]` |
| *L-Ergothioneine (mM)* 1 | 1.0 (0.1) `[M]` | 13 (0.0) `[M]` | |
| 5 | 1.2 (0.0) `[M]` | 9.3 (0.2) `[M]` | |
| 10 | 1.2 (0.1) `[M]` | 7.4 (0.1) `[M]` | |
| *Thiourea (mM)* 1 | 1.1 (0.1) `[M]` | 18 (0.1) `[M]` | |
| 5 | 1.1 (0.1) `[M]` | 15 (0.0) `[M]` | |
| 10 | 1.1 (0.2) `[M]` | 11 (0.3) `[M]` | |

ᵃ Fluorescence in relative units (±range, *n* = 2).
ᵇ AGE–BSA test was performed as indicated in Section 2.

*Note on the concentration units:* the **Rutin block is headed (µM), every other compound block is
headed (mM)**. The text layer collapses µ to m and would have read this as 200 mM; the page image
settles it.

### Numbers printed in figure captions and running text

| quantity | value | class | where |
|---|---|---|---|
| **k₂ₙd, D-penicillamine + phenylglyoxal** | **24.8 ± 1.3 /M s** | `[M]` | **Fig. 6 caption, p. 369** (printed twice, panels A and B) |
| **k₂ₙd, aminoguanidine + phenylglyoxal** | **0.4 ± 0.01 /M s** | `[M]` | **Fig. 6 caption, p. 369** |
| **k₂ₙd, L-cysteine + phenylglyoxal** | **0.63 ± 0.04 /M s** | `[M]` | **§3.4 running text, p. 368** |
| **k₂ₙd, D-penicillamine + phenylacetaldehyde** | **1.83 ± 0.05 /M s** | `[M]` | **Fig. 6 caption, p. 369** |
| D-penicillamine vs aminoguanidine, as the authors state it | "more than **60 times faster**" | `[M]` | §3.4, p. 368 |
| phenylglyoxal vs phenylacetaldehyde, as the authors state it | "approximately **14 times less efficiently**" | `[M]` | §3.4, p. 368 |
| AGE fluorescence signal, complete reaction vs H1 alone (1 mL geometry) | "about **20 times** greater than background" | `[M]` | §3.1, p. 365 |
| MALDI-TOF-MS of the methylglyoxal / D-penicillamine adduct | [M+Na]⁺ = **226 Da**, against a calculated **203.26 Da** for C₈H₁₃NO₃S | `[M]` / `[C]` | §2.6 p. 364 and §3.3 p. 368 |
| ¹H NMR, methylglyoxal adduct (D₂O, ppm) | 1.03 (3H, s, CH₃), 1.05 (3H, s, CH₃), 1.90 (3H, s, CO-CH₃), 3.96 and 3.98 (diastereomeric, 1H, s, CH-COOH) | `[M]` | §2.6, p. 364 |
| ¹³C NMR, methylglyoxal adduct | carbonyls at **δ ≈ 199 ppm (CH₃C=O)** and **δ = 209 ppm (COOH)** | `[M]` | §3.3, p. 368 |
| ¹H NMR, phenylglyoxal adduct (D₂O, ppm) | 1.38 (3H, s, CH₃), 1.45 (3H, s, CH₃), 4.12 (1H, s, CH-COOH), 7.42–7.82 (5H, m, ArH) | `[M]` | §2.7, p. 364 |
| conversion, phenylglyoxal + D-penicillamine, preparative | **>90 % in 40 min** at 10 mM / 20 mM, 50 mM KH₂PO₄ pH 7.4, room temperature | `[M]` | §2.7, p. 364 |
| screen result, best thiols | "even at the highest concentration tested (10 mM), L-Cys-OMe, L-Cys and NAC achieved only **50 % inhibition**" | `[M]` | §3.2, p. 367 |
| cell survival, fibroblasts, 1 mM D-penicillamine + 600 µM methylglyoxal | **47 % survival** | `[M]` | §3.5, p. 370 |
| cell survival, keratinocytes, same conditions | **100 % survival** | `[M]` | §3.5, p. 370 |
| growth inhibition threshold | **600 µM methylglyoxal completely growth inhibitory** in both cell lines | `[M]` | §3.5, p. 369 |

**Figure-only:** every point of Figs. 3, 4A, 6A, 6B, 7A–D and 8A–B, and all SDS-PAGE band intensities
in Fig. 4B. Per house rule, values read off those axes are not typed.

### Arithmetic on the printed constants (all mine)

**1. The scavenger ladder as within-study ratios.**

| pair | ratio (mine) | the paper's own words |
|---|---:|---|
| D-penicillamine / aminoguanidine on phenylglyoxal | 24.8 / 0.4 = **62.0×** | "more than 60 times faster" ✓ |
| D-penicillamine / L-cysteine on phenylglyoxal | 24.8 / 0.63 = **39.4×** | not stated |
| L-cysteine / aminoguanidine on phenylglyoxal | 0.63 / 0.4 = **1.58×** | "comparable" ✓ |
| phenylglyoxal / phenylacetaldehyde on D-penicillamine | 24.8 / 1.83 = **13.6×** | "approximately 14 times" ✓ |

Every one of my ratios reproduces the paper's own verbal statement, which is the only arithmetic
check this paper permits. **The α-oxo substitution is worth 13.6× on the same nucleophile** — that is
the paper's cleanest single-variable result, because the leaving structure, the buffer, the
temperature and the detection are identical between the two runs.

**2. What 0.63 M⁻¹ s⁻¹ means as a half-life in a pot (mine, and conditional).** At 37 °C in this
buffer, with cysteine at 1 mM, the pseudo-first-order constant on the α-dicarbonyl is
`0.63 × 1e-3 = 6.3e-4 s⁻¹`, i.e. a **half-life of about 18 minutes (mine)**. At 10 mM cysteine it is
about **110 s (mine)**. This is offered only to show that at physiological temperature the adduct
route is already fast on a Maillard timescale; **it is not a number to ship**, and it says nothing
about 140 °C (Flags 1).

**3. Table 1 read as a thiol ladder (mine).** Taking the complete reaction under air as 21 and each
thiol's 10 mM day-5 value: cysteamine 6.0 (**71 % suppression, mine**), L-Cys and L-Cys-OMe both 9.0
(**57 %**), NAC 9.4 (**55 %**), GSH 14 (**33 %**), thiourea 11 (**48 %**), ergothioneine 7.4
(**65 %**), D-penicillamine 1.4 (**93 %**), D,L-penicillamine 1.7 (**92 %**), aminoguanidine 1.8
(**91 %**). **GSH is the weakest thiol in the table**, and D,L-homocysteine (27) and L-Cys-Gly (40)
*raise* the signal above control, which the authors attribute to those compounds glycating on their
own account (§3.2, p. 367). These are fluorescence suppressions of a 5-day endpoint, **not rate
constants and not stoichiometries**, and they must not be converted into either.

## 4. Numbers the repository can use

Every row below shares: **10 mM phosphate buffer, pH 7.4, 37 °C, HPLC-followed disappearance of the
α-dicarbonyl at 254 nm, means ± SD of four measurements, atmosphere unstated** — except the Table 1
rows, which are the 5-day 50 mM KH₂PO₄ pH 7.4 37 °C histone-H1/ADP-ribose plate assay.

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **second-order constant, L-cysteine + phenylglyoxal** | **0.63 ± 0.04** | **M⁻¹ s⁻¹**, order 2 (1 in each reactant) | 37 °C, pH 7.4, 10 mM phosphate, 50 µM PG, 250/500 µM Cys | §3.4 running text, p. 368 | **`measured_rate`** |
| **second-order constant, D-penicillamine + phenylglyoxal** | **24.8 ± 1.3** | M⁻¹ s⁻¹, order 2 | as above, 250/500 µM D-Pen | Fig. 6 caption, p. 369 | **`measured_rate`** |
| **second-order constant, D-penicillamine + phenylacetaldehyde** | **1.83 ± 0.05** | M⁻¹ s⁻¹, order 2 | as above, with phenylacetaldehyde in place of PG | Fig. 6 caption, p. 369 | **`measured_rate`** |
| second-order constant, aminoguanidine + phenylglyoxal | **0.4 ± 0.01** | M⁻¹ s⁻¹, order 2 | as above | Fig. 6 caption, p. 369 | `measured_rate` (a **hydrazine**, not a thiol — carried as the paper's internal yardstick only) |
| activation energy for any of the above | **— none** | kJ/mol | single temperature (37 °C) only | — | **absent**; do not assign a barrier from this paper |
| α-oxo enhancement, same nucleophile | **13.6** | ratio | 37 °C, pH 7.4, D-penicillamine, phenylglyoxal vs phenylacetaldehyde | 24.8/1.83 (mine), matching the paper's "approximately 14 times" (§3.4, p. 368) | **`within_study_ratio`** |
| D-penicillamine / L-cysteine reactivity | **39.4** | ratio | 37 °C, pH 7.4, phenylglyoxal | 24.8/0.63 (mine) | `within_study_ratio` |
| D-penicillamine / aminoguanidine reactivity | **62.0** | ratio | as above | 24.8/0.4 (mine), matching "more than 60 times" | `within_study_ratio` |
| preparative conversion, PG + D-penicillamine | **>90 % in 40 min** | % | 10 mM PG, 20 mM D-Pen, 50 mM KH₂PO₄ pH 7.4, **room temperature** | §2.7, p. 364 | **`measured_bound`** (a lower bound on conversion; the temperature is not 37 °C) |
| AGE-fluorescence suppression at 10 mM thiol, 5 d | GSH 14, L-Cys 9.0, L-Cys-OMe 9.0, NAC 9.4, cysteamine 6.0, ergothioneine 7.4, thiourea 11, D-Pen 1.4, D,L-Pen 1.7, homocysteine 27, L-Cys-Gly 40, against a complete-reaction 21 (air) / 23 (argon) | relative fluorescence units | 1.5 mg/mL histone H1 + 1 mM ADP-ribose, 50 mM KH₂PO₄ pH 7.4, 37 °C, 5 days, 300 µL | Table 1, p. 366 | **`level_only`** — an endpoint fluorescence, not a rate, not a yield, not a stoichiometry |
| oxygen-independence of the glycation system | **23 (1.1) under argon + 5 mM DTPA vs 21 (0.9) under air**, day 5 | relative fluorescence units | as above | Table 1, p. 366 | **`within_study_ratio`** (**1.10× (mine)**) — the paper's demonstration that its α-dicarbonyl chemistry needs no O₂ |
| cell survival, 1 mM D-penicillamine against 600 µM methylglyoxal | **100 % (HaCat keratinocytes) / 47 % (CF3 fibroblasts)** | % of untreated growth | DMEM + 10 % FBS, 5 % CO₂, 37 °C, 72 h | §3.5, p. 370 | `level_only` — a cell-count endpoint |
| any rate, order or barrier for **thiol → disulfide** | **— none, anywhere in this paper** | — | — | — | **absent** |

### Adduct branch or redox branch?

**All four rate constants are ADDUCT branch. None is redox.** The assignment is not an inference from
stoichiometry — the paper isolates and structurally characterises the products:

| constant | branch | how it was decided |
|---|---|---|
| L-cysteine + phenylglyoxal, 0.63 M⁻¹ s⁻¹ | **ADDUCT** | Assigned by the paper's own framing: §3.4 is headed "Comparative reaction kinetics of α-dicarbonyl **trapping**", and cysteine is introduced as a member of the same thiazolidine-forming series ("the requirement of an α-amino-β-mercaptoethane structure", Discussion p. 370). L-cysteine has exactly the α-amino-β-mercapto motif that closes the ring. **Caveat, and it is real: the cysteine adduct itself is never isolated in this paper** — only the two D-penicillamine adducts are. The assignment is by homology within the paper, not by direct product characterisation for cysteine (Flags 6). |
| D-penicillamine + phenylglyoxal, 24.8 M⁻¹ s⁻¹ | **ADDUCT**, directly proven | Product isolated by preparative HPLC and assigned by ¹H and ¹³C NMR as **2-benzoyl-5,5-dimethyl-thiazolidine-4-carboxylic acid** (§2.7 p. 364, §3.3 p. 368): five aromatic protons, no aldehyde proton, the C-2 thiazolidine proton exchanged with D₂O. A C–S bond, not an S–S bond. |
| D-penicillamine + methylglyoxal (no rate constant; preparative only) | **ADDUCT**, directly proven | **2-acetyl-5,5-dimethyl-thiazolidine-4-carboxylic acid**, [M+Na]⁺ 226 Da against C₈H₁₃NO₃S calc. 203.26, ¹H NMR fully assigned, ¹³C carbonyls at 199 and 209 ppm, **no aldehyde proton** (§3.3, p. 368). |
| D-penicillamine + phenylacetaldehyde, 1.83 M⁻¹ s⁻¹ | **ADDUCT**, and explicitly the *reversible* kind | The paper states that monocarbonyl-derived thiazolidines "form reversibly with subsequent release of the aldehyde" (p. 371). This is an aldehyde, not an α-dicarbonyl, and it is present only as the α-oxo control. |
| aminoguanidine + phenylglyoxal, 0.4 M⁻¹ s⁻¹ | **ADDUCT**, and not a thiol at all | Aminoguanidine is a hydrazine; its product is a triazine/hydrazone, not a sulfur species. Carried only as the within-study yardstick. |

**The redox branch is absent by construction, not by oversight.** The paper's whole design goal is a
glycation assay that "proceeds independent of oxygen and therefore **excludes** identification of
inhibitory compounds acting as antioxidants" (Abstract, p. 361). The authors deliberately built a pot
in which redox chemistry cannot be the answer. **Anything a reader takes from this paper about
thiol oxidation would be taken from a system engineered to have none.**

**One thing this does say about the pre-registration's redox step, indirectly.** If the redox route
`α-dicarbonyl + 2 R-SH → hydroxyalkanone + RSSR` competed effectively with adduct formation at 37 °C
for cysteine, some fraction of the phenylglyoxal loss in Fig. 6 would have gone that way and the
isolated D-penicillamine adducts would not have been >90 % of the converted phenylglyoxal (§2.7).
**They were.** That is a *within-study bound* at 37 °C for D-penicillamine and, by the homology
argument only, a suggestion for cysteine. It is **not** a bound at 140 °C, where the redox route's
higher barrier (if it has one) could invert the ordering entirely (Flags 1).

## 5. Flags

1. **TEMPERATURE TRANSFER IS THE DOMINANT PROBLEM, AND THIS PAPER GIVES NO WAY TO PAY FOR IT.**
   Every constant here is at **37 °C**; the wave concerns **140 °C**, a gap of **103 K**. There is
   **only one temperature in the paper**, so no Arrhenius pair, no barrier, no way to extrapolate
   within the paper's own evidence. To transport 0.63 M⁻¹ s⁻¹ from 310 K to 413 K one must import an
   activation energy from somewhere else, and the size of the resulting move is enormous and almost
   entirely determined by the imported number: at Ea = 50 kJ/mol the factor is **≈ 1.1e2**; at
   80 kJ/mol **≈ 1.5e3**; at 122.2 kJ/mol (`ZHANG_EA_THIOL_TO_DISULFIDE_KJ_MOL`, the barrier the
   repository already carries on `ch_dimer_*`) **≈ 1.1e5** (all **mine**, Arrhenius,
   `exp[−Ea/R(1/413 − 1/310)]`). **Three orders of magnitude of the answer are chosen by the
   modeller, not measured by Wondrak.** That is the whole cost, stated plainly: this paper fixes a
   *number* and leaves its *temperature dependence* entirely free, which is the same failure mode
   the pre-registration identifies in Whitfield. Any use of these constants above ~40 °C is a
   `derived_assumption`, never a `measured_rate`.
2. **The probe is phenylglyoxal, and the pot's α-dicarbonyls are not.** Phenylglyoxal was chosen for
   its UV chromophore, not for its relevance. It is an **aryl** α-oxoaldehyde and it is the most
   electrophilic of the family; the repository's α-dicarbonyls (2,3-pentanedione, 2,4-pentanedione,
   3,4-hexanedione in the norfuraneol pot; deoxypentosones on the sugar lane) are **alkyl diketones
   and deoxyosones with no aldehyde carbon at all**. The paper's own α-oxo control shows that
   removing one oxo group costs **13.6×**; nothing here says what replacing an aryl aldehyde with a
   symmetric alkyl diketone costs, and a 2,3-diketone has no aldehydic carbon for the Schiff base
   that Fig. 5 requires. **The cysteine constant should not be attached to 2,3-pentanedione without
   an explicit, labelled assumption.**
3. **Methylglyoxal and glyoxal have no rate constant here.** They appear only preparatively (24 h,
   0.20 M buffer) and in cell culture (72 h). Anyone reaching for "Wondrak's methylglyoxal constant"
   will not find one.
4. **The two reactant ratios are pooled, so the pseudo-first-order assumption is not independently
   auditable.** §2.8 says the 5:1 and 10:1 values "were in good agreement" and prints only the
   pooled result. The demonstration that k₁ₛₜ scales with scavenger concentration is explicitly
   "(data not shown)". The **order in the scavenger is asserted, not shown**, in the published
   record.
5. **Irreversibility is asserted on a literature citation and a stability observation, not measured.**
   §3.3 is headed "Irreversible trapping" and the Discussion says the α-dicarbonyl adduct "is stable
   in water" (ref. [44]). **No hydrolysis rate, no equilibrium constant and no back-reaction is
   measured for the α-dicarbonyl thiazolidines.** If the repository models the adduct as an
   irreversible sink on this paper's authority, that authority is a citation.
6. **The cysteine constant's branch assignment is by homology, not by product isolation.** Only the
   D-penicillamine adducts were characterised. Cysteine could in principle also form the open-chain
   hemithioacetal without ring closure, and the two would be indistinguishable in an assay that
   follows only α-dicarbonyl disappearance. **What is certain is that the 0.63 M⁻¹ s⁻¹ is not a
   disulfide-forming rate**, because no disulfide was formed, sought or reported anywhere.
7. **Nothing here is at a Maillard pH.** pH 7.4 throughout (the preparative phenylglyoxal run is at
   pH 7.4 too). The Whitfield pot the pre-registration must reproduce is at **pH 4.5**. Thiol
   nucleophilicity is thiolate-gated: at pH 7.4 cysteine (thiol pKa ≈ 8.3) is a few per cent
   ionised, at pH 4.5 it is ~1e-4 of that. `ph_state.py` already carries a thiolate fraction, so the
   repository can in principle pay this, but **the correction is roughly four orders of magnitude
   and it must be applied explicitly**, not absorbed.
8. **The buffer is 10 mM phosphate for the kinetics and 50 mM or 0.20 M for the preparative and
   screening pots.** Phosphate catalyses enolisation and dicarbonyl chemistry. The Whitfield pot is
   **0.5 M phosphate** (per the pre-registration's own charge correction), a 50× higher phosphate
   than the kinetic runs here. No buffer-concentration series is reported.
9. **The kinetic runs' atmosphere is not stated.** The argon/DTPA control belongs to the *screening*
   assay only. If any of the phenylglyoxal loss in Fig. 6 were oxygen-assisted, nothing in §2.8
   would have caught it.
10. **The rate constants live in a figure caption.** There is no rate-constant table, no supporting
    information on disk, and no raw kinetic data. The four numbers are exactly as printed and cannot
    be cross-checked against anything.
11. **Table 1 has n = 2 and reports a range, not a standard deviation.** It is a screen. Treating any
    Table 1 value as a quantitative measurement of thiol consumption would be a category error: it
    is the fluorescence of a *protein* after five days, and a compound can lower it by trapping
    dicarbonyls, by chelating, by quenching (which the AGE–BSA column is there to catch) or by
    glycating in the compound's own right (which is what L-Cys-Gly and homocysteine visibly do).
12. **GSH's weakness in Table 1 is not a rate statement.** GSH at 10 mM leaves 14 of 21 — the poorest
    suppression of any simple thiol here. It is consistent with Zheng 2022's finding that GSH's
    reaction with α-dicarbonyls is limited under physiological conditions, and it is **the same
    class of evidence** (an endpoint level, not a rate), so the two corroborate each other only
    weakly.
13. **Conflict of interest is declared and is material to the compound ranking.** Two authors are
    principals in Niadyne Inc., which sponsored the research; the paper's conclusion is that the
    D-penicillamine pharmacophore should be developed as a therapeutic. This does not touch the
    arithmetic, but the *choice* of which compounds got kinetic characterisation (D-penicillamine,
    and two comparators) rather than a full ladder is a commercial choice, and it is why GSH, NAC,
    cysteamine and ergothioneine have no constants.
14. **What to request from the authors:** (i) the k₁ₛₜ values at 250 and 500 µM separately, for both
    ratios and all compounds; (ii) any temperature series at all — even one extra temperature would
    convert every constant here from a point to a barrier; (iii) product characterisation for the
    **L-cysteine** + phenylglyoxal adduct; (iv) whether any disulfide (cystine, penicillamine
    disulfide) was ever detected in the HPLC traces and simply not reported; (v) constants for
    methylglyoxal and glyoxal rather than phenylglyoxal.
15. **What this paper does NOT contain:** any disulfide; any GSSG; any temperature other than 37 °C
    (kinetics) or room temperature (one preparative run); any pH other than 7.4; any activation
    energy; any alkyl α-diketone; any constant for GSH; any volatile or headspace measurement; any
    food system; any DFT or other computation (so nothing was excluded under the no-DFT policy).

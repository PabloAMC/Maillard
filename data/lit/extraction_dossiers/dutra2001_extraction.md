# Dutra, Knudsen, Curi & Bechara 2001 — EXTRACTION
DOI: 10.1021/tx015526r (checked against Crossref, 2026-11-01)
### *"Aerobic Oxidation of Aminoacetone, a Threonine Catabolite: Iron Catalysis and Coupled Iron Release from Ferritin"*
### Aminoacetone (AA, the α-aminoketone that the engine's B18 Strecker step makes from methylglyoxal + glycine) autoxidises in Chelex-treated 0.1 M phosphate pH 7.4, **37 °C**, to **methylglyoxal + NH4⁺ in 1 : 1 stoichiometry with the O2 consumed**. Oxygraph rates: apparent k2 = 0.160 M⁻¹ s⁻¹; Fe(II)EDTA ×3, Cu(II) aquo ×20; SOD −90 %; phosphate 50 -> 100 mM ×2.2; pH 5.8 -> 8.6 ×~80. EPR spin trapping. **37 °C only; no product formation under N2 / argon (the N2 arm is only for ferritin iron release); no other buffer.**

**Read 2026-10-31 from `data/articles/dutra2001.pdf`** (7 pp., Chem. Res. Toxicol. 14:1323-1329; clean text layer).
Figs. 1, 2, 4 and Tables 1-2 checked against 200 dpi renders of pp. 1325-1326 (read error ~ ±0.5 µM/min on Fig. 4,
~ ±5e-5 s⁻¹ on Fig. 2). Figs. 5-7 (EPR, ferritin, heme proteins) read for content only.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Questions: triage item 7 (oxygen as a reagent) and item 5 (pyrazines; the aminoketone pool), plan list item 9
(phosphate). The engine's pyrazine lane (`src/kinetic_core/network.py`, B18 / B22) makes aminoacetone from MGO + amino
acid and condenses it to 2,5-dimethylpyrazine; it has **no oxidative aminoacetone sink**.

| need | what Dutra 2001 gives | tag |
|---|---|---|
| O2 / argon contrast | O2 is the measured reagent (Clark electrode). Product: MGO **205 ± 13 µM** and NH4⁺ **212 ± 15 µM** from 5 mM AA in air-equilibrated buffer (~200 µM O2), 2 h — "stoichiometric" with O2 consumed [M]. No N2 / argon product arm; N2 purging blocks ferritin iron release [M] | [M] / [!] |
| phosphate vs other buffer | **None at one pH.** Phosphate 50 vs 100 mM: **1.14 vs 2.47 µM O2/min** (5 mM AA, pH 7.4; "not shown") [M]. The pH profile (Fig. 4) runs in 100 mM phosphate, so pH co-varies with the HPO4²⁻ fraction [!] | [M] / [!] |
| metals / chelators | Table 2 (5 mM AA, 37 °C, µM O2/min): none 3.9; Fe(II)EDTA 30 µM 12.5; Fe(III)EDTA 30 µM 3.9; Fe(III)EDTA 200 µM 15.5; Fe(II)ATP 150 µM 5.0; Fe(III)ATP 3.8; Fe(II)citrate 300 µM 5.3; Fe(III)citrate 6.1. Cu(II) aquo 50 µM ×20; Cu(II)EDTA 100 µM no effect. DTPA and Desferal up to 3 mM: no effect on the uncatalysed rate (no adventitious metal after Chelex) [M] | [M] |
| rate constants | **k2 = 0.160 ± 0.007 M⁻¹ s⁻¹** (slope of kobs vs [AA], 3-20 mM, kobs = rate / [O2], [O2] assumed 200 µM) [M]. Fit has a positive x-intercept (~3 mM) [Z] | [M] / [Z] |
| Ea | none (37 °C; EPR at 25 °C) | — |
| isotope labels | none | — |
| transfer to 100-140 °C | **Not shown.** Physiological chemistry; the authors make no food or temperature claim | [!] |

**Answer.** For the model this is the **measured rate law and product of the aerobic aminoketone -> dicarbonyl
back-reaction** (AA + O2 -> MGO + NH4⁺ + H2O2), superoxide-chain propagated, accelerated by Fe(II) and Cu(II) and by
phosphate, steeply pH dependent (deprotonated AA, pKa 6.8). In a closed food pot the sink is capped by dissolved O2
(~0.2 mM at 37 °C in air, less hot), so it can matter only where O2 is resupplied or aminoketones are sub-mM; it then
competes with the O2-requiring dihydropyrazine -> pyrazine oxidation [Z]. Confidence that the reaction exists in
aerobic Maillard pots at pH ≥ 7: ~0.8; that the 37 °C k2 transfers to 100-140 °C within ×10: low (~0.3; no Ea, chain
kinetics, metal-dependent) [Z]. Use: a declared rival / band for an "aminoketone + O2" edge, not a constant.

## 2. Conditions

| item | as printed |
|---|---|
| aminoacetone | hydrochloride (Hepworth), recrystallised; stock 500 mM in deaerated water; pKa of the α-ammonium **6.8** (NaOH titration, Fig. 4 inset) |
| buffer | 100 mM phosphate pH 7.4, Chelex-100 treated (all buffers); 50 mM in one comparison; heme-protein work in 10 mM (text) / 50 mM (Fig. 7 caption) Tris-glycine pH 7.8 [!] |
| temperature | 37 °C (oxygraph, ferritin); 25 °C (EPR); 36 °C (Fig. 7) |
| O2 | YSI 53 oxygraph, Clark electrode, air-equilibrated (~200 µM assumed); MGO product run "in a flask with minimal headspace", 2 h |
| analysis | MGO: 1,2-benzenediamine derivatisation (60 °C, 3 h) after HClO4 stop, HPLC-DAD 315 nm; NH4⁺: enzymatic kit; EPR with DMPO 100 mM |
| replicates | triplicates; rates "reproducible within 15 %" |

## 3. Quantitative results

**Kinetics** [M]: k2 = 0.160 ± 0.007 M⁻¹ s⁻¹. Fig. 2 [F]: kobs ≈ 1.3e-4 / 3.0e-4 / 6.2e-4 / 1.0e-3 / 2.7e-3 s⁻¹ at
3 / 5 / 7 / 10 / 20 mM AA; the line crosses zero near 3 mM, and the 3-mM point lies above it [Z]. Fig. 1: sigmoidal O2
uptake (lag ~3 min, autoxidation-chain shape); 5 mM AA consumes the dissolved O2 in ~50 min, with Fe(II)EDTA in ~18 min
[F].

**Scavengers, Table 1** (5 mM AA, % inhibition of O2 uptake) [M]: catalase 53 (+Fe(II)) / 40 (−Fe); SOD 22 (+Fe) /
90 (−Fe); semicarbazide 80 (±Fe). Xanthine / xanthine oxidase (O2•⁻ source) accelerates O2 uptake by 1 mM AA (Fig. 3).

**Phosphate** [M]: 1.14 (50 mM) vs 2.47 (100 mM) µM O2/min, 5 mM AA, pH 7.4: apparent order ~1.1 in phosphate [Z].

**pH, Fig. 4** (5 mM AA + 30 µM Fe(II)EDTA, 100 mM phosphate, 37 °C; µM O2/min) [F]: ~0.3 (5.8), 0.5 (6.0), 0.7 (6.2),
1.1 (6.4), 1.8 (6.6), 4.0 (6.8), 7.7 (7.0), 9.0 (7.2), 16 (7.4), 20 (7.6), 20.5 (7.8), 21.7 (8.0), 21 (8.2), 22.3 (8.4),
24 (8.6). Inflection 7.2 [M]. Between pH 5.8 and 8.6 both the free-base AA fraction (pKa 6.8) and the HPO4²⁻ fraction
(pKa2 ~6.8 at 37 °C, I ~0.2) rise ×~11 [Z], so the profile cannot separate substrate deprotonation from dianion
catalysis.

**Derived** [Z]: with O2 held at air saturation, AA's pseudo-first-order loss would be k2·[O2] ≈ 3.2e-5 s⁻¹ (t½ ≈ 6 h)
at 37 °C, pH 7.4, 100 mM phosphate, uncatalysed; ×3 with 30 µM Fe(II)EDTA. In a sealed pot, total conversion ≤ the
O2 inventory (1 O2 per AA).

## 4. Mechanism evidence and its kind (Bechara laboratory, São Paulo)

- **Stoichiometry / mass balance (kind 5):** MGO : NH4⁺ : O2 ≈ 1 : 1 : 1 [M].
- **Trapping (kind 3):** EPR spin trapping (DMPO-•OH from DMPO-OOH; a 6-line carbon-centred adduct assigned to the AA
  enoyl radical [P assignment, "?" in the text]); DMSO -> DMPO-•CH3 confirms •OH with Fe(II)EDTA.
- **Within-study factors (kind 6):** [AA], metal complex and redox state, SOD / catalase / semicarbazide, phosphate
  50 vs 100 mM, pH 5.8-8.6.
- **Proposed [P]:** phosphate-catalysed enolisation of deprotonated AA (by analogy with 5-aminolevulinic acid, ref 12),
  one-electron reduction of O2 by the enol, iminoacetone hydrolysis to MGO + NH4⁺; with Fe(II), a concerted
  two-electron transfer to H2O2 then Fenton (Scheme 3).

## 5. Bearing on open questions

| item | verdict | why |
|---|---|---|
| Triage 7, oxygen as a reagent | **SUPPORTS** a new O2 edge: aminoketone + O2 -> dicarbonyl + NH4⁺ | measured stoichiometry and rate, 37 °C |
| Triage 5, pyrazines (aminoketone partition) | **adds a rival sink** to "aminoketone -> pyrazine vs -> pyrrole" (evidence map 2.12 "pyrazines: pyrrole sink SINGLE") | aerobic loss back to MGO; its weight at 100-140 °C in a sealed vial is O2-limited [Z] |
| Phosphate (triage 9) | **SUPPORTS** phosphate acceleration of an α-carbonyl enolisation-gated step (×2.2 for ×2 phosphate) | one pH; no other buffer; pH profile confounded with HPO4²⁻ fraction |
| Metals as a hidden factor (thornalley1999, smuda2010 / voigt2010 DTPA confound) | **SUPPORTS** that Fe / Cu at tens of µM multiply α-carbonyl autoxidation ×3-20, while chelated Fe(III) and Cu(II)EDTA do not | chelator identity matters: EDTA keeps Fe(II) active, DTPA / Desferal do not add to the uncatalysed rate |
| CF-011 / CF-076, hub split | **LEAVES UNCHANGED** | no sugar chemistry |
| CF-073 (Strecker yield from MGO + Phe at pH 3) | **LEAVES UNCHANGED** | AA oxidation is negligible at pH 3 by Fig. 4's trend (≤ 0.3 µM/min already at 5.8) [Z] |

## 6. Disagreements found

1. Internal [!]: the uncatalysed rate at 5 mM AA, 100 mM phosphate, pH 7.4, 37 °C is **3.9 µM O2/min** in Table 2 but
   **2.47 µM O2/min** in the phosphate comparison (text, "not shown"); same nominal conditions, 37 % apart (stated
   reproducibility 15 %). Fig. 2's kobs at 5 mM (3.0e-4 s⁻¹ × 200 µM = 3.6 µM/min) matches Table 2.
2. Internal [!]: the abstract says AA 1-7 mM; the k2 text says 3-20 mM (Fig. 2 spans 3-20). The Tris-glycine buffer is
   10 mM in the Methods and 50 mM in the Fig. 7 caption.
3. No other dossier on disk measures aminoacetone + O2; no cross-lab check possible (single lab).

## 7. Ledger-ready measurements (two readers)

1. k2 = 0.160 M⁻¹ s⁻¹ with the Fig. 2 points (5) — the only measured aminoketone autoxidation rate on disk; flag the
   x-intercept and 37 °C.
2. MGO / NH4⁺ / O2 stoichiometry (printed) — the product branch of the edge (MGO, not acetol or pyruvate).
3. Table 2 metal-complex rates (8 rows) and Table 1 scavenger inhibitions — metal and chain factors.
4. Fig. 4 pH profile (15 points) — keyed with "phosphate dianion co-varies".
5. Phosphate 50 / 100 mM pair — keyed with the Table 2 inconsistency.

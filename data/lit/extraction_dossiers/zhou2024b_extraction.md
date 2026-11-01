# Zhou, Huang, Cui, Chen, Hayat, Zhang & Ho 2024 (JAFC 72:5878) — EXTRACTION
### *"Exogenous Alanine Promoting the Reaction between Amadori Compound and Deoxyxylosone and Inhibiting the Formation of 2-Furfural during Thermal Treatment"*
### Purified Ala-xylose Amadori compound (AX-ARP, 20 mM) +- L-alanine 20 mM (14N or 15N), unbuffered water, initial pH 7.0 (falls during heating; Fig. S2, not on disk), 120 °C, 30-120 min: ARP (HPLC-ELSD), 3-DX / 1-DX / GO / MGO (OPD, ex situ), furfural and 2,3-pentanedione (HS-SPME), 14N- and 15N-alanine (UPLC-MS/MS), and an ARP + 3-DX condensation product (m/z 336, peak area only). **ONE initial pH; no temperature series; the condensation product is not quantified.**

**Read 2026-10-31 from `data/articles/zhou2024b.pdf`** (9 pp., born-digital ACS, clean text layer; stamped "Downloaded
... on 03 September 2026"). Prose from `pdftotext -layout`; Figs. 1, 2, 5 and 6 read from 150 / 200 dpi page renders
(read error ~ +-0.03 mM on Fig. 2a, +-0.01 mM on Figs. 2b-d, +-0.3 mM on Figs. 1b and 5, +-0.002 mM on Fig. 1a). Figs. 3-4
(mass spectrum, scheme) checked on the render. The Supporting Information (Figs. S1-S4: chromatograms, **pH time
course S2**, condensation-product peak areas S3, MGO + Ala pot S4) is NOT on disk.

**Naming.** This is the paper `zhou2024_extraction.md` calls "Zhou 2024a (JAFC 72(11), 5878)" -- the method reference
of that dossier's pH-8 AX-ARP +- 15N-Ala arm. `zhou2024` (10.1021/acs.jafc.4c03706, JAFC 72:18630) is a different
paper; this file is NOT a duplicate.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Tong Zhou, Meigui Huang, Heping Cui, Pusen Chen, Khizar Hayat, Xiaoming Zhang* (Jiangnan University, Wuxi), Chi-Tang Ho* (Rutgers) |
| Venue | *J. Agric. Food Chem.* **2024**, 72 (11), 5878-5886 |
| DOI | **10.1021/acs.jafc.4c00021** (printed on every page; Crossref 2026-10-31: same title, seven authors, vol. 72, issue 11, pp. 5878-5886, issued 2024-03-11) |
| Dates | received 2 Jan 2024; revised 22 Feb 2024; accepted 29 Feb 2024; published 11 Mar 2024 |
| Laboratory | **Jiangnan University (Zhang X / Cui H) with Ho (Rutgers): the "Ho network" of the evidence map** -- same laboratory as Zhou 2022 / 2023 / 2023c / 2024 / 2025b, Zhai 2019-2023, Wei 2022, Cui 2022, Hou 2026 |

## 1. What this paper answers for the model

| model need | what Zhou 2024b gives | tag |
|---|---|---|
| (a) 1-DX / 3-DX split as a FLUX or a standing level | **Standing levels only** (OPD added after the quench). 3-DX / 1-DX (mM / mM) = **2.6 / 5.6 / 7.8 / 9.7** at 30 / 60 / 90 / 120 min without Ala, **2.2 / 5.2 / 7.9 / 9.6** with Ala: the ratio RISES ~4x in time while the pH falls (Fig. S2, not on disk) | [F] / [Z] |
| (b) deoxyosone lifetimes | Not measured. Bounds from the 90-120 min decline against the ARP loss (assumes 3-DX / 1-DX made only from the ARP): **3-DX ~50 min to ~5 h, 1-DX ~7 to ~100 min at 120 °C** (pH drifting below 7) | [Z] |
| (c) amine catalysis of enolisation / cleavage | Exogenous Ala speeds ARP loss (apparent first-order k 0-30 min: **0.0139 vs 0.0190 /min, x1.37**; 0-60 min x1.21) but does NOT raise any deoxyosone or GO / MGO, and LOWERS the release of the ARP's own alanine per ARP consumed (**~0.69 vs ~0.91-1.0 mol/mol**). No sign that the extra ARP loss runs through enolisation | [F] / [Z] |
| (d) a base-favoured 3-DX sink | A **3-DX sink proportional to [ARP]**: an ARP + deoxyxylosone condensation product (m/z 336.1306, [M+H]+ of ARP + DX - H2O) detected in both pots, more with Ala (Fig. S3, peak area only). Ala lowers 3-DX by 12-19 % [F]. **No pH series**, so no pH law; no cysteine | [M] / [F] / [P] |
| (e) observation operator | OPD **ex situ** (0.4 mL sample + 0.4 mL reagent, 4 h dark); 3,4-hexanedione internal standard; **no authentic 3-DX / 1-DX standard** -> semi-quantitative (CF-062 class). GO / MGO from purchased 40 % standards | [M] |

**Answers, weighed.**
1. **Flux vs level.** Nothing here is a formation split. The standing ratio climbs from ~2.5 to ~10 in 90 min in both
   pots, as in every other ARP pot of this laboratory (Hou 2026 2.6 -> 10.9 with the pH measured; Zhou 2022 5 -> 22;
   Wei 2022 1.3 -> 2.2). The flat-in-time ratio of signature (i) is a fed-TTCA observation (Zhai 2019), not a property of
   ARP pots [Z]. Confidence the climb is real and not an artefact of the semi-quantitative scale: ~0.8 (a shared
   response-factor error would shift the ratio, not its slope).
2. **Alanine is not an enolisation catalyst here** (my reading, ~0.7). The extra ARP loss with Ala (~1.9 mM at 30 min)
   is not matched by more 3-DX, 1-DX, GO or MGO, and 14N-Ala release per ARP consumed falls by ~30 % [Z]. The authors'
   mechanism -- ARP condenses with 3-DX, Ala accelerates it -- is MS identification plus a peak area [P]; the mass of
   the product is equally that of an ARP + 1-DX adduct (same formula), which the authors acknowledge by taking 3-DX
   "as an example" [Q][!].
3. **A 3-DX sink, but not base-tested.** The condensation is a second-order 3-DX exit (k [ARP][3-DX]) that the engine
   does not have. Its pH law is untested here; Zhou 2024 (pH 8, same pot design) is the only same-laboratory second pH
   point and its Fig. 3 is figure-only in that dossier.
4. **The authors' "equilibrium shifted toward ARP" story contradicts their own Fig. 1b** [!]: they write ARP
   degradation "was suppressed to some extent" [V], while the ARP falls faster with Ala at every time.

## 2. Conditions [M]

| item | as printed |
|---|---|
| ARP synthesis | Ala : Xyl 1 : 2 in water, pH 7.5, 80 °C 60 min, vacuum dehydration 80 °C 20 min, Dowex 50WX8 (H+), water then 0.18 M ammonia; ELSD-pure fraction freeze-dried; identity by UPLC-ESI-MS and NMR (as ref 19). **Purity not stated** [!] |
| Model pots | ARP 20 mmol/L; ARP + 14N-Ala 20 mmol/L; ARP + 15N-Ala 20 mmol/L. **No buffer** (water); initial pH 7.0 with 2 M NaOH; **not held**; pH measured after reaction (Fig. S2 -- "pH ... was <7", not on disk) |
| Vessel | pressure-resistant bottles (100 mL); fill volume, headspace gas, deaeration **not stated** |
| Heating | oil bath **120 °C**, 30 / 60 / 90 / 120 min; ice-water quench |
| Side pot | Ala 20 + MGO 20 mmol/L, pH 7.0, 120 °C, 0-120 min (Fig. S4, not on disk) |
| ARP | HPLC-ELSD, XBridge BEH amide, external standard (purified ARP); RT 21.3 min |
| Dicarbonyls | **ex situ**: 0.4 mL reaction solution + 0.4 mL reagent "1% OPD with DTPA dissolved in HEPES (1 mmol/L)" [V; 1 mM HEPES is likely a typo for 1 mol/L, cf. Zhou 2022 "1 M HEPES"] + 10 uL 3,4-hexanedione (86 mmol/L) internal standard; 4 h dark; HPLC-DAD 315 nm, SunFire C18; RT 3-DX 15.6, 1-DX 16.1, GO 18.2, MGO 24.1 min |
| 14N / 15N-Ala | UPLC-TQD, MRM, external standards: 14N y = 3.25e6 x + 31080 (R2 0.9907), 15N y = 4.53e6 x + 22207 (R2 0.9985), x in mmol/L |
| Volatiles | HS-SPME DVB/CAR/PDMS, 3 mL + 2 mL sat. NaCl, IS 1,2-dichlorobenzene, 60 °C 30 min; DB-Wax; furfural y = 0.0191 x - 0.3498 (R2 0.9900), 2,3-pentanedione y = 0.0139 x + 0.0669 (R2 0.9985), x in ug/L |
| Replicates | n = 3; ANOVA / Duncan |

## 3. Quantitative results

**Fig. 1b, ARP (mmol/L), 120 °C [F]:**

| min | 0 | 30 | 60 | 90 | 120 |
|---|---|---|---|---|---|
| ARP | 20 | 13.2 | 8.0 | 5.2 | 4.3 |
| ARP + Ala | 20 | 11.3 | 6.6 | 4.3 | 3.9 |

**Fig. 2, mmol/L [F]:**

| min | 3-DX | 3-DX +Ala | 1-DX | 1-DX +Ala | GO | GO +Ala | MGO | MGO +Ala |
|---|---|---|---|---|---|---|---|---|
| 30 | 1.38 | 1.21 | 0.54 | 0.545 | 0.12 | 0.135 | 0.07 | 0.065 |
| 60 | 2.10 | 1.78 | 0.375 | 0.34 | 0.185 | 0.185 | 0.12 | 0.10 |
| 90 | 2.02 | 1.74 | 0.26 | 0.22 | 0.31 | 0.35 | 0.185 | 0.185 |
| 120 | 1.85 | 1.49 | 0.19 | 0.155 | 0.41 | 0.455 | 0.255 | 0.24 |

**Fig. 1a, furfural (mmol/L) [F]:** ARP 0.013 / 0.030 / 0.0425 / 0.045; +Ala 0.0085 / 0.024 / 0.0345 / 0.0375 at 30-120
min. Furfural at 120 min = **0.29 % of ARP consumed** (0.23 % with Ala) [Z]: furfural is a minor 3-DX exit here.

**Fig. 5, alanine (mmol/L) [F]:**

| min | 14N-Ala (ARP pot) | 14N-Ala (ARP + 15N-Ala) | 15N-Ala (ARP + 15N-Ala) |
|---|---|---|---|
| 0 | 0 | 0 | 20 |
| 30 | 7.8 | 6.0 | 17.5 |
| 60 | 12.2 | 9.2 | 17.2 |
| 90 | 13.8 | 10.9 | 17.1 |
| 120 | 14.3 | 11.3 | 17.0 |

[Z] 14N-Ala released per ARP consumed: ARP pot **1.15 / 1.02 / 0.93 / 0.91** (30-120 min; the 30-min value exceeds
1 by more than the read error [!]); with 15N-Ala **0.69 / 0.69 / 0.69 / 0.70**. Exogenous 15N-Ala lost: 2.5-3.0 mM.

**Fig. 6 [F]:** A420 ARP 0.265 / 0.33 / 0.40 / 0.445, +Ala 0.33 / 0.45 / 0.54 / 0.595; 2,3-pentanedione (ug/L) ARP
33.2 / 16.8 / 6.6 / 3.8, +Ala 43 / 31.5 / 11.5 / 7.8 (30-120 min).

**Fig. 3a [M]:** condensation product m/z 336.1306 with fragments 318.1194, 300.1092 (-H2O, -2H2O), 222.1009 (ARP+H),
204.0874 (base peak), 90.0567 (Ala+H). Peak-area time course only in Fig. S3 (not on disk): "reached a maximum amount
at 30 min" [V], stable to 150 min [Q]. **Pyrazines: none detected** in either pot [Q].

**Derived lifetimes [Z]** (90-120 min window, both pots similar; ARP loss 0.030 mM/min, 3-DX -0.006 mM/min at ~1.9 mM,
1-DX -0.0023 mM/min at ~0.23 mM): k3 between 0.0032 /min (formation = 0) and 0.019 /min (all ARP loss -> 3-DX);
k1 between 0.010 and 0.14 /min. In-situ pH in that window unknown (< 7).

## 4. Mechanism evidence and its kind

| claim | kind | strength |
|---|---|---|
| ARP + deoxyxylosone condensation product exists | K3-ish (MS identification of a product, no standard, no quantity) | weak-moderate: exact mass fits ARP + C5H8O4 - H2O; 3-DX vs 1-DX not distinguished |
| Ala accelerates the condensation | peak area (Fig. S3, not on disk) | weak: one-pot comparison, no response factor needed (same analyte), but not seen by me |
| Exogenous Ala does not condense with ARP or DX | negative MS search | moderate (no LOD) |
| Endogenous Ala release is lowered by exogenous Ala | **K1 isotope labelling (15N-Ala)** with calibrated MRM | **strong** for the observation; the interpretation (equilibrium shift) is [P] |
| 3-DX not -> furfural, not -> GO / MGO with Ala | within-study factor (+-Ala) | moderate |
| 3-DX / 1-DX standing levels | measured intermediate, ex situ OPD, semi-quantitative | moderate for shape, weak for scale |

**Laboratory:** Jiangnan University (Zhang X, Cui H) with Ho (Rutgers).

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| Evidence map 2.1, "amine ACCELERATION of 2,3-enolisation / beta-cleavage: NONE" | **LEAVES UNCHANGED** (adds a mild negative) | Ala speeds ARP loss but no deoxyosone or C2 / C3 fragment rises; 1-DX is unchanged. Single laboratory |
| Evidence map 2.2 / 2.12, the 3-DX sink | **SUPPORTS a non-thiol 3-DX sink exists** (ARP + 3-DX), single lab, MS + peak area | a second-order sink k [ARP][3-DX] would shorten 3-DX lifetime early (ARP high) and is not in the engine |
| Prereg sec. 19-22 (v4 feasible "base-favoured 3-DX sink"; v5 Cys trap rejected) | **LEAVES UNCHANGED on the pH law** | no pH series; the condensation's pH dependence is unknown. Zhou 2024 (pH 8) Fig. 3 could give a two-point law |
| Hub H1 (both osones quasi-steady; engine 3-DX lives ~35 h) | **SUPPORTS** | 3-DX lifetime <= ~5 h, plausibly ~1 h, at 120 °C [Z] |
| Hub H1's "ratio flat in time" premise / signature (i) | **CONTRADICTS for ARP pots** | ratio rises ~4x in 90 min here (and in three other same-lab ARP pots); signature (i) is a TTCA-pot observation |
| Hub H7 (pH drift in unbuffered pots) | **SUPPORTS (qualitatively)** | pH falls below 7 (Fig. S2, not on disk); Hou 2026 measures the same drift in a GX-ARP pot |
| CF-060 (Jiangnan AX-ARP dicarbonyl scale spread) | **EXTENDS** | a fifth AX-ARP-alone 120 °C pot, see sec. 6 |
| CF-062 (deoxyosones semi-quantitative) | **SUPPORTS** | internal standard, no authentic 3-DX / 1-DX standard |
| CF-071 (Cys double ratio, Zhai 2019 vs 2023) | LEAVES UNCHANGED | no cysteine. Alanine's double ratio (3-DX/1-DX with Ala over without) is 0.87 / 0.94 / 1.02 / 0.99 [Z]: an ordinary amino acid does NOT reproduce signature (iii) |

## 6. Disagreements found

1. **Same laboratory, same nominal pot, 3x apart (extends CF-060).** AX-ARP 20 mM, water, initial pH 7, 120 °C,
   pressure bottle, 30-120 min, same ARP synthesis and OPD method (ref 17 = Zhou 2022):

   | quantity | Zhou 2022 (zhou2022_extraction.md Fig. 1) | Zhou 2024b | ratio |
   |---|---|---|---|
   | 3-DX (mM) | 4.3 / 6.1 / 5.9 / 5.8 | 1.38 / 2.10 / 2.02 / 1.85 | ~0.33 |
   | 1-DX (mM) | 0.82 / 0.61 / 0.33 / 0.26 | 0.54 / 0.375 / 0.26 / 0.19 | ~0.7 |
   | GO (mM) | 0.47 / 0.75 / 0.83 / 1.12 | 0.12 / 0.185 / 0.31 / 0.41 | ~0.3 |
   | MGO (mM) | 0.10 / 0.14 / 0.18 / 0.27 | 0.07 / 0.12 / 0.185 / 0.255 | ~0.9 |
   | released Ala (mM) | 2.2 / 4.2 / 5.8 / 5.0 | 7.8 / 12.2 / 13.8 / 14.3 | **~2.9** |
   | furfural (umol/L) | 4.6 / 8.6 / 10.9 / 13.1 (443-1261 ug/L) | 13 / 30 / 42.5 / 45 | **~3** |
   | ARP (mM) | 12.0 / 8.3 / 6.9 / 6.2 | 13.2 / 8.0 / 5.2 / 4.3 | ~1 early, 0.7 late |

   Calibrated species (Ala by different instruments: HPLC-ELSD in 2022, UPLC-MS/MS in 2024b; furfural by SPME external
   curves in both; GO by purchased standard) differ by ~3x too, in opposite directions, so a deoxyosone response factor
   alone (CF-062) cannot explain it. Candidate causes: ARP batch (purity unstated in both); amine quantification method
   (ELSD vs MRM: the 2024b release exceeds stoichiometry at 30 min, 2022's stays at ~36 %); headspace / O2 (unstated in
   both); a 4-h vs 4-h OPD protocol is the same, so not the derivatisation time.
2. **Released Ala exceeds ARP consumed at 30 min (1.15 mol/mol)** [Z][!]: internal inconsistency at the level of
   +-0.3 mM read error on both curves; an ARP under-recovery or Ala over-recovery of ~15 %.
3. **Text vs figure on ARP stability** (sec. 1, item 4) [!].
4. **Engine:** no ARP + 3-DX condensation or any [ARP]-proportional 3-DX exit (`src/kinetic_core`: TDP exits are
   `r_tdp_fur`, `r_tdp_mgo` (off), `r_tdp_cys_trap` (off)).

## 7. Ledger-ready measurements (two readers)

1. Fig. 5 14N / 15N-Ala (calibrated MRM, isotope-resolved): endogenous amine release per ARP consumed, +-exogenous Ala
   -- a mass-balance row for the ARP -> DX + amine step, response-factor-immune in ratio.
2. Fig. 2a / 2b 3-DX and 1-DX +-Ala (semi_quantitative, shape only): the ratio's time slope in an ARP pot and the
   +-Ala double ratio.
3. Fig. 1b ARP +-Ala (calibrated): ARP lifetime and the amine factor on it at pH 7 initial.
4. Fig. 2c / 2d GO / MGO (calibrated per CF-062's adjudication of this method) and Fig. 1a furfural (calibrated).
5. Not keyable until fetched: Fig. S2 (the pH time course, both pots) -- the in-situ-pH input every row above needs.

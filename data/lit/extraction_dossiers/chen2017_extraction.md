# Chen et al. 2017 — EXTRACTION (free oleic acid heated at 140 C with 0 / 2 / 5 % water; DMPO spin-trapping EPR; HS-SPME-GC-MS/MS volatiles as area %)
### The free-oleic-acid volatile slate at 140 C, with 2-decenal, 2-undecenal and nonanal assigned to the 8- and 9-hydroperoxides, and a drawn 2-decenal -> 2,4-decadienal step.

**Source on disk:** `data/articles/chen2017.pdf` (owner's download, 2026-09-08; "Article in press",
Food Chemistry 2017). Read from the `pypdf` text layer; Table 1 (page 6) is clean and was checked
against a rendering of the page. Scheme 1 (page 5) and Fig. 4 (page 6) are drawn mechanisms,
described in words. Figs 1, 3, 5 are plots: FIGURE-ONLY (the plateau numbers quoted in the text are
recorded).

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of water content on thermal oxidation of oleic acid investigated by combination of EPR spectroscopy and SPME-GC-MS/MS" |
| Authors | Hongjian Chen, Peirang Cao, Bo Li, Dewei Sun, Yong Wang, Jinwei Li, Yuanfa Liu (Jiangnan Univ.; Jinan Univ.) |
| Venue | Food Chemistry 2017 (PII S0308814616318313) |
| DOI | 10.1016/j.foodchem.2016.11.008 |
| Naming | "8-COOH" and "9-COOH" in Section 3.4 are the paper's (mis)notation for the oleate 8- and 9-hydroperoxides. "RO" = rapeseed oil (used only for AV/PV, Fig 1); "WC" = water content. |
| Cited for the drawn step | Warner, Neff, Byrdwell & Gardner 2001, JAFC 49:899 (2-decenal -> 2,4-decadienal in heated triolein) |

## 1. Why it matters

Free oleic acid (not an oil) heated at a frying-adjacent temperature with the whole volatile slate
printed: it gives the oleate aldehyde set (nonanal, octanal, decanal, 2-decenal, 2-undecenal) as
relative percentages, plus the alcohols (1-octanol, 1-heptanol, pentanol) and acids that the R-radical
partners of each scission become, plus 2-pentylfuran. It states the 8-OOH -> 2-undecenal and 9-OOH ->
2-decenal / nonanal / octanoic acid assignments and draws a second-generation step, 2-decenal ->
2,4-decadienal, that no other dossier carries. The water-content effect (EPR radical yield up 140 % at
5 % water; 2,4-decadienal share doubles) is a process-variable handle. It is the only paper in this
batch with a stated argument that free acid is more oxidisable than the triglyceride.

## 2. Methods as they matter to a model

- **Substrate:** oleic acid from J & K Chemical (Shanghai). Purity: the sentence reads "Oleic acid and
  DMPO (99 % purity)"; whether 99 % applies to the oleic acid is ambiguous. The volatile slate (hexanal
  2.6 %, 2-heptenal 3.1 %, 2,4-decadienal 2.8 %) indicates a linoleate impurity.
- **Volatiles run:** 2.00 mL oleic acid in a sealed 20 mL headspace vial, oil bath **140 C / 30 min**;
  water added to the vial to 0, 2, 5 % (w/w presumably; not stated). HS-SPME DVB/CAR/PDMS 50/30 µm,
  equilibrate >= 5 min, extract **50 C / 40 min** with magnetic stirring. GC-MS/MS TSQ Quantum XLS,
  HP-5MS 30 m x 0.25 mm x 0.25 µm, He 0.8 mL/min, 50 C (3 min) -> 125 C at 8 C/min -> 165 C at 4 C/min
  (3 min) -> 230 C at 10 C/min (2 min); EI 70 eV, scan 50-400. Identification: RI vs pure references
  and NIST/Wiley match.
- **Quantification: area percent of total volatiles** (Table 1 values are "%"; footnote "-" =
  undetected). No internal standard, no calibration. Class totals for 0 % water: aldehydes 49.85 %,
  "others (acids, furanic compounds etc.)" 37.56 %, hydrocarbons 3.70 %, alcohols 8.02 %, ketones 1.36 %
  (data not shown); aldehydes 59.97 % at 2 % water, 56.78 % at 5 % water.
- **EPR:** 100 µL oleic acid + 20 µL DMPO in toluene (150 mM; final ~25 mM), 4 mm quartz tube, cavity
  at **140 C**, spectra every 2 min to ~15 min, dark; Bruker EMXplus 9.85 GHz, 20 mW, 1.0 G modulation.
  Spin adduct assignments: alkyl adduct aN = 14.89 G, aHb = 20.74 G, g = 2.00669 (dominant); alkoxyl
  adduct aN = 13.64 G, aHb = 7.73 G, aHc = 1.50 G, g = 2.00672 (early, then destroyed); "PBN-like"
  adduct aN = 14.88 G, aH = 1.81 G, g = 2.00685 (constant). Signal intensity is the "DApp" parameter of
  the carbon-centred adduct, arbitrary units. n = 3.
- **Oleic acid remaining:** FAME/GC-FID, ratio to C13:0 internal standard ("Au", Fig 5, FIGURE-ONLY;
  text: degradation rose steeply 0 -> 2 % water, slightly 3 -> 5 %).
- **Rapeseed oil (76.5 % oleate, alumina-purified, tocopherol-free):** AV and PV at 140 C vs water
  content (Fig 1, FIGURE-ONLY); text: PV 2.680 -> 3.875 mg/kg (sic, unit as printed) from 1 % to 3 %
  water (+44.59 %).
- **Cited kinetics:** Chen, Lee & Schanus 1992: methyl linoleate hydroperoxide decomposition
  decelerated by 0 -> 2 % water (hydrogen bonding) — the opposite direction to this paper's finding at
  140 C.

## 3. Tables re-typed

### Table 1. "the major volatile compounds of fried oleic acids at 140 C with different moisture contents." Values are area % of total volatiles at 0 / 2 / 5 % water; "-" = undetected.

| volatile | 0 % | 2 % | 5 % |
|---|---:|---:|---:|
| heptane | 0.63 | 0.62 | 0.51 |
| octane | 2.41 | 1.65 | 1.16 |
| pentanal | 0.61 | 0.67 | 0.42 |
| hexanal | 2.58 | 2.93 | 1.55 |
| heptanal | 1.30 | 1.61 | 1.03 |
| **octanal** | 2.34 | 2.73 | 1.96 |
| **nonanal** | 9.08 | 11.27 | 7.50 |
| **decanal** | 0.60 | 0.90 | 1.06 |
| (E)-2-hexenal | 0.18 | 0.18 | 0.12 |
| (E)-2-heptenal | 3.11 | 4.08 | 2.34 |
| (E)-2-octenal | 1.03 | - | - |
| 2-nonenal | 1.65 | 2.80 | 2.07 |
| **(E)-2-decenal** | 11.94 | 14.44 | 14.92 |
| **(E)-2-undecenal** | 12.12 | 12.82 | 16.57 |
| **(E,E)-2,4-decadienal** | 2.78 | 3.63 | 5.78 |
| pentanol | 0.76 | 0.87 | 0.54 |
| 1-heptanol | 1.75 | 2.04 | 1.24 |
| 1-octanol | 3.03 | 3.06 | 2.75 |
| formic acid | 0.98 | 0.71 | 0.55 |
| hexanoic acid | 0.73 | 0.69 | 0.88 |
| heptanoic acid | 2.70 | 1.76 | 2.51 |
| octanoic acid | 8.81 | 4.43 | 7.07 |
| nonanoic acid | 1.35 | 0.65 | 1.05 |
| n-decanoic acid | 11.32 | 7.05 | 8.54 |
| dodecanoic acid | 1.68 | 1.08 | 1.63 |
| **2-pentylfuran** | 1.40 | 2.21 | 1.23 |
| 5-ethyldihydro-2(3H)-furanone | 0.24 | 0.24 | 0.24 |

Column sums of the listed rows: 87.1 / 84.9 / 84.9 %; the remainder is unlisted minor peaks.

Text numbers on Table 1: "(E)-2-Decenal (11.94 %), (E)-2-Undecenal (12.12 %) and nonanal (9.08 %)" the
most abundant aldehydes at 0 % water; 2-undecenal at 5 % water "36.72 % higher than that of the sample
without water" (16.57/12.12 = 1.367); 2-decenal "increased 20.94 % when the water content rose from 0 %
to 2 %" (14.44/11.94 = 1.209); 2,4-decadienal at 5 % water "nearly twice" the 0 % value (5.78/2.78 =
2.08).

EPR numbers in the text (Fig 3, arbitrary intensity): 2 % water samples 0.90-1.05 from 2 to 12 min vs
0.55-0.75 without water; 5 % water rose to a maximum 2.20 at ~10 min then plateaued at 1.80 after 12
min; abstract: "EPR intensity plateau of the samples with 5 % water content was 140 % higher than the
samples without water".

### Scheme 1 — "Structure formulas of formation of lipid-derived radical adducts" (drawn; described)

Oleic-type chain R-CH=CH-CH2-CH2-R' + O2 (catalyst) -> allylic peroxyl R-CH=CH-CH(OO•)-CH2-R' (LOO•);
LOO• + DMPO -> DMPO-OOL adduct; the DMPO-OOL adduct decomposes to a short-chain alkoxyl R''-CH=CH-CH(O•)-R'''
(RO•) + DMPO-oxide; RO• + DMPO -> DMPO-OR adduct; alkyl radical R-CH=CH-CH•-CH2-R' + DMPO -> DMPO-L
adduct. The scheme says nothing about C-C scission; it is the spin-trap bookkeeping.

### Fig. 4 — "The proposed formation mechanism of (E,E)-2,4-decadienal" (drawn; described)

(E)-2-decenal (CHO-CH=CH-CH2-(CH2)5-CH3) -> loses H• from C4 (the allylic methylene) -> C4 radical
-> + O2, H• -> 4-hydroperoxy-2-decenal -> "- H2O2" -> (E,E)-2,4-decadienal (CHO-CH=CH-CH=CH-(CH2)4-CH3).
Text: "(E)-2-Decenal may undergo a hydroperoxidation at 2-decenal allylic methylene carbon by classical
free radical mechanism to form a hydroperoxide group. A double bonds would be produced by further loss
of hydrogen peroxide to produce (E,E)-2,4-decadienal (Fig. 4) (Warner, Neff, Byrdwell, & Gardner,
2001)." The step is a net H2O2 elimination as drawn; the paper does not give a radical-level
mechanism for the elimination.

## 4. Routes and numbers the repository can use

All levels: area % of total headspace volatiles, oleic acid, 140 C / 30 min, sealed 20 mL vial, SPME
50 C / 40 min; replicate count for Table 1 not stated (EPR and FAME n = 3).

| route | reactant -> product | mechanism as drawn / stated | measured numbers with units and conditions | evidence class |
|---|---|---|---|---|
| OL-8-A | oleate 8-OOH -> (E)-2-undecenal | stated: "(E)-2-Undecenal was formed by the beta-scission of 8-COOH" (no figure) | 12.12 / 12.82 / 16.57 % at 0 / 2 / 5 % water | measured_level (area %); mechanism stated, not drawn |
| OL-9-A | oleate 9-OOH -> (E)-2-decenal + octanoic acid | stated: "(E)-2-Decenal, nonanal and octanoic acid ... were generated by the breakdown of 9-COOH. Therefore ... 9-COOH was the main hydroperoxide" | 2-decenal 11.94 / 14.44 / 14.92 %; octanoic acid 8.81 / 4.43 / 7.07 % | measured_level; mechanism stated |
| OL-9-B | oleate 9-OOH -> nonanal | stated (same sentence) | nonanal 9.08 / 11.27 / 7.50 % | measured_level; mechanism stated |
| OL-11-B, OL-8-B | -> octanal, decanal | not discussed; present | octanal 2.34 / 2.73 / 1.96 %; decanal 0.60 / 0.90 / 1.06 % | measured_level |
| OL-ALC | alkyl-radical partners -> 1-octanol, 1-heptanol, octane, heptane | not drawn; consistent with 10-OOH / 11-OOH A-scission partners (Cao 2020 nomenclature) | 1-octanol 3.03 / 3.06 / 2.75 %; 1-heptanol 1.75 / 2.04 / 1.24 %; octane 2.41 / 1.65 / 1.16 %; heptane 0.63 / 0.62 / 0.51 % | measured_level |
| OL-DEC-DIEN | **(E)-2-decenal -> (E,E)-2,4-decadienal** | allylic C4 H-abstraction, O2, 4-OOH, - H2O2 (Fig 4) | 2,4-decadienal 2.78 / 3.63 / 5.78 %; ratio 2,4-decadienal / 2-decenal 0.23 / 0.25 / 0.39 | mechanism_drawn (second-hand, Warner 2001); measured_ratio (within-study, water dependence) |
| OL-PF | oleic acid (with linoleate impurity) -> 2-pentylfuran | not discussed | 1.40 / 2.21 / 1.23 % | measured_level; source ambiguous (see flag 1) |
| OL-H2O | water 0 -> 5 % raises radical yield and shifts slate to 2-alkenals / dienal | EPR + Table 1 | EPR plateau +140 % (5 % vs 0 %); aldehyde class 49.85 -> 59.97 -> 56.78 % | measured_ratio (arbitrary units) |

## 5. Rule sketches (repository suggestions, not the paper's)

Keys: `nonanal`, `hexanal`, `2_pentylfuran`, `heptanal`, `e_2_octenal` exist; **octanal, decanal,
(E)-2-decenal, (E)-2-undecenal, (E,E)-2,4-decadienal, 1-octanol, 1-heptanol** do not
(`DECADIENAL` exists in `structures.yml` as `CCCCC/C=C/C=C/C=O`). No oleic acid or oleate hydroperoxide
isomer structures exist (only the lump `LOOH_OL`); see cao2020_extraction.md §5 for proposed isomer
SMILES and for the oleate scission rules S1/S2, which this paper's assignments (8-OOH -> 2-undecenal,
9-OOH -> 2-decenal / nonanal / octanoic acid) corroborate on the free acid.

**S1. 2-alkenal -> 2,4-alkadienal (new, net; second-generation oxidation of an aldehyde).**
Pattern: `O=[CH1:1][CH1:2]=[CH1:3][CH2:4][CH2:5][#6:6]` -> `O=[CH1:1][CH1:2]=[CH1:3][CH1:4]=[CH1:5][#6:6]`
(net loss of H2 as H2O2 with O2). Applies to any 2-alkenal with a CH2-CH2 beyond C3.
- positive control: `CCCCCCC/C=C/C=O` ((E)-2-decenal) -> `CCCCC/C=C/C=C/C=O` (DECADIENAL, (E,E)-2,4-decadienal)
- second positive: `CCCCCCCC/C=C/C=O` ((E)-2-undecenal) -> `CCCCCC/C=C/C=C/C=O` (2,4-undecadienal)
- negative control: nonanal `CCCCCCCCC=O` (saturated aldehyde: no C2=C3, must not fire); 2-pentenal `CC/C=C/C=O` (no C5 methylene beyond C4... has C4 CH2 and C5 CH3 only — a rule needing CH2-CH2 must not fire) ; hexanal.
- caveat: the paper's own data cannot separate this route from direct linoleate-derived 2,4-decadienal (linoleate impurity present); the evidence for the step is Warner 2001 (not on disk).

**S2. Oleate 8-OOH / 9-OOH A-scission (cao2020 S2)** — this paper is a second, free-acid witness for
the 2-undecenal / 2-decenal + octanoic acid products; positive controls as in cao2020 S2 with the free
acid: `CCCCCCC/C=C/C(OO)CCCCCCCC(=O)O` -> `CCCCCCC/C=C/C=O` + `CCCCCCCC(=O)O`.

**S3. Alkyl radical -> alcohol / alkane (R• + •OH -> ROH; R• + R'H -> RH).** Not drawn here; the
1-octanol / octane and 1-heptanol / heptane pairs are the measured witnesses. If a net rule is written
for the 10-OOH / 11-OOH A-scission partner, the volatile products should be written as the alcohol and
the alkane, not as an aldehyde. Aldehyde -> alcohol (e.g. hexanal -> 1-hexanol) is NOT supported by
this paper (1-hexanol not detected).

## 6. Flags

1. **Substrate purity.** Hexanal, 2-heptenal, 2-octenal, 2-nonenal, (E,E)-2,4-decadienal, pentanal and
   2-pentylfuran in the slate are the linoleate signature; the paper attributes 2,4-decadienal to
   2-decenal oxidation without a linoleate-free control. Any within-study ratio involving these
   compounds carries the impurity.
2. **Area %, no internal standard, no response factors**; percentages of a total whose denominator
   changes with water content. Only within-column rank order and cautious cross-column ratios are
   usable.
3. **Table 1 replicate count is not stated.** EPR and FAME were "in triplicate"; the SPME table has no
   SD.
4. **"8-COOH / 9-COOH"** notation and "beta-scission of 8-COOH" for 2-undecenal: with Cao 2020's
   geometry this is the C7-C8 cleavage (Cao's "A-scission"); the two papers use different words for the
   same bond. The dossier records the bond, not the label.
5. **Fig 4's "- H2O2" elimination** from a 4-hydroperoxy-2-alkenal is a net arrow; a 4-hydroperoxide
   would more usually go to a 4-oxo- or 4-hydroxy-2-alkenal. Treat as mechanism_drawn (second-hand).
6. **EPR intensities are arbitrary units** of a spin-adduct steady state at 140 C in the presence of
   25 mM DMPO in toluene; "140 % higher" is a plateau ratio, not a rate.
7. **The sealed 20 mL vial with 2 mL acid** has ~18 mL air (~0.16 mmol O2 vs ~6.3 mmol oleic acid):
   oxygen-limited after the early phase (the paper says so for the EPR tube: "the system was like under
   oxygen-depleted condition").
8. **PV unit printed as "mg/kg"** (should be meq/kg or mmol/kg); values 2.680-3.875 kept as printed.
9. No DFT, no rate constants; nothing to mark inadmissible.

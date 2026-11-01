# Weenen & Tjan 1992 — EXTRACTION
DOI: 10.1021/bk-1992-0490.ch017 (checked against Crossref, 2026-11-01)
### *"Analysis, Structure, and Reactivity of 3-Deoxyglucosone"* (ACS Symp. Ser. 490, *Flavor Precursors*, ch. 17, pp. 217-231)
### 3-DG purified on an α-cyclodextrin HPLC column; in water it is ~5 major + >= 6 minor isomers with NO detectable C=O (13C-NMR): almost entirely (hemi)acetal / hydrate. Stable in water (70 °C, 2 h: no change); 2 N AcOH 100 °C 1 h -> 40-50 % HMF; base -> metasaccharinic acid. Methylpyrazines from fructose / glucose / 3-DG + asparagine in 1,2-propanediol + K2HPO4, 130 °C, 1 h; [1-13C]glucose labelling of the pyrazine methyls (~26 %) and ring CH (~20 %). **No Strecker aldehyde is measured, and there is NO cysteamine trapping in this chapter.**

**Read 2026-10-31 from `data/articles/weenen1992.pdf`** (15 pp., ABBYY OCR of a scan, ACS Publications download;
DOI printed in the page margin). Text from `pdftotext`, every page checked against 110 dpi renders (scratchpad
`w92-*.png`; Table 1 and the experimental subscripts, garbled in the OCR, read from the renders). Figs. 1-3
(HPLC traces, UV spectrum, 13C-NMR spectra) hold no quantities used here; Schemes 1-8 are mechanisms.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

The brief and `docs/guides/EVIDENCE_MAP_2026-10-31.md` sec. 4 (line ~799) list this chapter as "cysteamine trapping
of 3-DG (K3, 2.2)", taken from `weenen1998_extraction.md` line 114, which cites "Weenen & Tjan **1994** with
cysteamine trapping; NMR". **This 1992 chapter contains no cysteamine, no trapping, and no Strecker-aldehyde
measurement** (confidence ~0.98, every page read). The cited trapping work is a different (1994) publication, not
on disk under this stem.

| model need | what Weenen & Tjan 1992 gives | tag |
|---|---|---|
| 3-DG as a Strecker donor (yield, rate) | **none** | [M] (absence) |
| why 3-DG is a weak Strecker donor vs MGO (CF-073's donor ratio) | **structural evidence:** no C=O carbon detectable by 13C-NMR; IR shows "only a very small absorption" at 1700-1800 cm^-1; ~5 major + >= 6 minor isomers (cyclic hemiacetals / hemiketals, bicyclic forms proposed); conversion to the free C-2 carbonyl "very slow" | [M] / [Q] / [P] (the bicyclic structures) |
| 3-DG lifetime in water | no change by 1H / 13C-NMR after 2 h at 70 °C in D2O; after months only "small amounts of HMF" | [Q] |
| 3-DG in acid | 1 N CF3COOD, 60 °C, 1 h: < 10 % to HMF; 2 N AcOH, 100 °C, 1 h: **40-50 % to HMF**; 1 N CF3COOD / D2O 100 °C 1 h: complete decomposition to HMF + formic acid; D at C-3 of HMF only 20-30 % | [M] |
| 3-DG in base | fast to metasaccharinic acid (0.1 M NaOH or 0.02 M Ca(OH)2); 0.1 M K2DPO4 / D2O pH 8: H-D exchange at C-3 complete in ~4 h at ambient; 0.1 M K2HPO4 in 1,2-propanediol (1 % water), 130 °C, 1 h: all 3-DG gone, **35-40 % metasaccharinic acid** | [M] |
| pyrazine carbon origin | [1-13C]glucose + Asn: ~26 % (20-31 %) 13CH3 in dimethylpyrazine; ~20 % (14-26 %) 13CH in all pyrazines; read as MGO from 1-DG vs from 3-DG | [M] / [P] (the assignment) |
| O2 and the dihydropyrazine oxidation | presence / absence of O2 did not change yield or distribution of methylpyrazines (fructose + Asn) | [Q] |

**Answers.**
1. **For CF-073** (MGO / 3-DG Strecker donor ratio 35 at Quest vs ~1.6 at DFA): the chapter, from the same Quest
   group as Weenen 2001, shows 3-DG is almost entirely in non-carbonyl cyclic / hydrated forms and opens slowly. That
   is a plausible reason for 3-DG's low Strecker yield in BOTH labs, not a reason the two labs disagree.
   Re-reading CF-073's numbers this way [Z]: the 3-DG arm AGREES across labs (Weenen 2001 PA 0.6 % at pH 3.2 vs
   Hofmann pH 3 PA ~0.7 %, from the prereg's PA + PAA 0.92 mmol/L at PA : PAA ~3 : 1; ~0.06 dex), while the MGO arm
   disagrees (21 % vs ~1.1 %, 1.28 dex). The conflict is in the C2 / C3 donors (MGO, GO), not the C6 osone.
2. **For the hub (3-DG lifetime, 2.1 / 2.2):** neutral water is slow (no NMR change at 70 °C / 2 h); acid and base
   are fast. Qualitative and partly in non-aqueous media; no rate constant.
3. **For pyrazines (evidence map 2.7, carbon origin):** a Quest K1 result (13C) at 130 °C in propanediol with
   Asn, adding a third lab beside Nestlé (Amrani-Hemaimi) and Reading (Low), with a different medium.

## 2. Conditions [M]

- **HPLC:** Cyclobond-III (α-cyclodextrin on silica), acetonitrile / water 88 / 12; RI (LOD 0.15 mg) or UV 200 nm
  (LOD 0.01 mg); preparative purification of > 100 mg. HMF on a BioRad fruit-quality column (0.002 N H2SO4);
  metasaccharinic acid by anion chromatography with pulsed electrochemical detection (0.1 N NaOH), standard
  synthesised from 3-methoxyglucose.
- **NMR:** 1H / 13C, DEPT, COSY, J-resolved at 100 / 200 / 360 MHz in CDCl3 or D2O; 13C/12C ratios in duplicate from
  13C satellites.
- **Pyrazines:** hexose + asparagine + K2HPO4 in 1,2-propanediol, 130 °C, 1 h; Likens-Nickerson extraction; GC with
  methoxypyrazine IS. 2,5- and 2,6-DMP not resolved ("dimethylpyrazine"). Concentrations not printed.
- **[1-13C]3-DG** from [1-13C]glucose by the Khadem method; **3-DG** by the Khadem / Madson-Feather bis(benzoylhydrazone)
  route.
- Replicates: duplicate 13C ratios; otherwise not stated.

## 3. Quantitative results

**Table 1 — "Formation of pyrazines from hexoses and asparagine"** [M] (units not printed; footnote: yields "not
optimized"; read as % [!]):

| C-source | N-source | pyrazine | methylpyrazine | dimethylpyrazine | trimethylpyrazine |
|---|---|---|---|---|---|
| fructose | Asn | - | 2.3 | 4.9 | 0.4 |
| glucose | Asn | - | 1.8 | 2.6 | 0.1 |
| 3-DG | Asn | - | 0.9 | 2.8 | 0.1 |

Ratios [Z]: DMP / MMP 2.1 (fructose), 1.4 (glucose), 3.1 (3-DG). "Very small amounts (< 0.05 %)" of 2,3-DMP and
tetramethylpyrazine [M].

**Other printed numbers** [M]: 3-DG -> HMF < 10 % (1 N CF3COOD, 60 °C, 1 h); 40-50 % (2 N AcOH, 100 °C, 1 h); D
incorporation at HMF C-3 20-30 % (vs 75-87 % from ARPs in 2 N AcOH / D2O, Feather & Russell 1969 [C]);
metasaccharinic acid 35-40 % (130 °C, 1 h, 0.1 M K2HPO4 / propanediol); 13CH3 in DMP ~26 % (20-31 %); 13CH in all
pyrazines ~20 % (14-26 %); HPLC LOD 0.15 mg (RI), 0.01 mg (UV 200 nm); UV maximum 302 nm (very small).

## 4. Mechanism evidence and its kind

Laboratory: **Quest International, Bussum (Weenen group)** — the same group as weenen1998 / weenen2001.
- **Isotope labelling (kind 1):** [1-13C]glucose -> pyrazines; [1-13C]3-DG -> HMF + formic acid (formic acid C and
  HMF aldehyde C both from C-1 of 3-DG) [M].
- **Directly measured structure (kind 4):** 13C-NMR / IR of purified 3-DG: no free carbonyl [M].
- **Within-study factor (kind 6):** acid vs neutral vs base on 3-DG fate; O2 presence vs absence on pyrazines [Q].
- **Proposed [P]:** bicyclic hemiacetal / hemiketal structures 11-14; HMF via a furylium ion without the free C-2
  carbonyl (Scheme 3); retro-aldol fragmentation of 1-DG / 3-DG to pyruvaldehyde, glyceraldehyde etc. (Schemes 5-8);
  dihydropyrazine oxidation by a dicarbonyl (MGO -> acetol) rather than O2 (Scheme 4).

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| CF-073 (Weenen 2001 vs Hofmann 2000, x14-20) | **LEAVES UNCHANGED, but localises it** | gives a structural reason for low 3-DG reactivity in both labs; re-read with CF-073's own numbers, the 3-DG yields agree (~0.06 dex) and the disagreement is MGO / GO only [Z]. Candidate causes then point to C2 / C3 donor handling (commercial aqueous MGO / GO speciation at 26 vs 100 mmol/L; prereg sec. 10.6) or stripping, not to a donor-generic quantification bias |
| Forced-input layer: kappa_TDG declared, 3-DG as donor | **CONTEXT** | slow ring opening supports a 3-DG Strecker flux limited by its open-chain fraction; no number |
| Forced-input layer SINGLE-LAB | **LEAVES UNCHANGED** | no Strecker measurement |
| Hub 3-DG lifetime (2.1 / 2.2) | **CONTEXT** | water at 70 °C stable for 2 h; consistent with ghiron1988's 100 °C losses within an Arrhenius factor [Z] |
| Pyrazine carbon origin (2.7, "CAN") | **SUPPORTS** (third lab, different medium) | K1 13C result, propanediol / 130 °C |
| Pyrazine oxidation step and O2 (2.8 context) | **weak SUPPORT** for an O2-independent dihydropyrazine oxidation | stated, data not shown |
| CF-074 / CF-075 | **LEAVES UNCHANGED** | Hofmann-internal |

## 6. Disagreements found

1. **`docs/guides/EVIDENCE_MAP_2026-10-31.md` (sec. 4, Weenen & Tjan 1992 row: "Cysteamine trapping of 3-DG (K3,
   2.2)") and the dossier brief** vs the PDF: no cysteamine trapping anywhere in this chapter. The source of the claim,
   `weenen1998_extraction.md:114`, cites "Weenen & Tjan **1994**". The K3 attribution to this stem is wrong; the 1994
   work is not on disk under `weenen1992`. (Hidalgo & Zamora 2004 is the on-disk paper that uses cysteamine trapping,
   for PA.)
2. Not a disagreement, a consistency: ghiron1988 (3-DG alone, water, 100 °C: ~33 % lost in 1 h after a fast initial
   drop) vs "no significant change" at 70 °C / 2 h here. With Ea ~100-130 kJ/mol the 100 -> 70 °C factor is ~15-30
   [Z], i.e. ~3-6 % expected at 70 °C / 2 h: within NMR's "no significant change". Ghiron's fast first-15-min drop
   (~30 %) is not mirrored by any instability here, which supports reading it as an equilibration / work-up effect.

## 7. Ledger-ready measurements

- **Table 1** (methylpyrazines from fructose / glucose / 3-DG + Asn, 130 °C, propanediol): low value (units unstated,
  non-aqueous); a shape comparator for the pyrazine sibling only.
- **13C incorporation** (26 % 13CH3 in DMP; 20 % ring 13CH) from [1-13C]glucose: a K1 carbon-origin number for the
  pyrazine lane's MGO-source split (1-DG vs 3-DG), medium caveat.
- **3-DG -> HMF 40-50 % (2 N AcOH, 100 °C, 1 h)** and **-> metasaccharinic acid 35-40 % (130 °C, phosphate /
  propanediol)**: fate shares for a hub check, both outside the engine's aqueous pH window.
- Nothing for the Strecker layer.

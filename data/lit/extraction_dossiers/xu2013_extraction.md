# Xu et al. 2013 (FRI) — EXTRACTION (xylose 1 : 1 [¹³C₅]-labelled, 26 mmol/L total, + L-cysteine 19 mmol/L + thiamine ± one beef-derived octa-/nona-/decapeptide, water, pH 4.8, 140 °C, 90 min, stirred Parr reactor; END-POINT volatiles by purge-and-trap GC-O-MS as relative concentrations; CAMOLA isotopomer splits; non-volatiles by ESI-Q-TOF)

### No TTCA, no Amadori compound, no cysteine and no time axis anywhere in this paper. It is in the brief because it shares the cysteine-xylose system; its only sulfur-lane content is a CAMOLA claim that MFT and FFT carry **no xylose label** — which contradicts Cerny & Briffod 2007 and Zhao 2019 for FFT and should not be used.

**Source on disk:** `data/articles/xu2013.pdf` (8 pp., typeset Elsevier version). Prose read from
`pdftotext`; **Tables 1, 2 and 3 re-read from 150 dpi page renders** (the text layer drops the
column structure of Table 2; scratchpad `img/xu13p-2.png`, `xu13p-3.png`, `xu13p-4.png`). Figs. 1-3
are a pathway scheme, an EI mass spectrum with a pathway, and ESI-Q-TOF spectra — not digitised (no
kinetic content). No Supporting Information is referenced.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation mechanism of volatile and non-volatile compounds in peptide–xylose Maillard reaction" |
| Authors | Qianqian Xu, Jianbin Liu, **Huanlu Song** (corresponding, songhl@th.btbu.edu.cn), Tingting Zou, Ye Liu, Songpei Zhang |
| Affiliation | Laboratory of Molecular Sensory Science, Beijing Technology and Business University (BTBU) — the Song group, not the Xie group of Cao 2017 / Zhao 2019 |
| Venue | Food Research International **54 (2013) 683-690**. Received 30 Apr 2013, accepted 28 Jul 2013 |
| DOI as printed | **10.1016/j.foodres.2013.07.066** ✔ matches the brief |
| Funding | NSFC 31171645 |
| Protocol parent | Lv & Song 2010 (not on disk) |
| Cited by | Zhao 2019 (`zhao2019_extraction.md`) for furfural from 3-deoxypentosone and 2-methylfuran from furfural reduction |

## 1. Why it matters (and why it barely does)

The paper's goal is peptide-specific Maillard products (cyclic dipeptides, a peptide-xylose
cross-link at m/z 320) from three synthetic beef peptides. For this engine:

1. **Nothing on TTCA, Amadori compounds or cysteine kinetics.** Cysteine is detected in the LC-MS
   ([M+H]⁺ 122.0689 "existed in all systems") but not quantified; no intermediate is followed; one
   time point.
2. **A carbon-origin claim that conflicts with the corpus.** "2-furfurylthiol and 2-methyl-3-furanthiol
   did not appear labeled molecules, and this was different from the former reports ... they may be
   generated from thiamine". For **MFT** that is directionally consistent with Cerny & Briffod 2007
   (`cerny2007_extraction.md`: MFT 75-90 % unlabelled in [¹³C₅]-xylose + Cys + thiamine — i.e.
   mostly thiamine-derived) and with the engine, which has a thiamine route (`r_thi_hmp` → HMP →
   `r_hmp_mft`). For **FFT** it **contradicts** Cerny & Briffod 2007 (FFT **> 99 % ¹³C₅**, i.e.
   exclusively xylose-derived, at pH 4-6) and Zhao 2019 (FFT 72 % M+5 from [¹³C₅]-TTCA), and it
   contradicts the engine's topology (FFT forms only from `FUR` + H2S; furfural in this very paper
   is 56 : 44 unlabelled : M+5, i.e. xylose-derived). The MFT/FFT signals are the smallest in Table 2
   (MFT 0.042-0.243), so the absence of a labelled cluster is most plausibly a detection-limit
   artefact. **Do not encode.**
3. **Pyrazines absent at pH 4.8 / 140 °C** ("pyrazines ... had not been found in such low pH and
   temperature conditions") — a qualitative negative consistent with the engine's base-catalysed
   pyrazine chemistry, not a number.

### How it relates to the companion dossiers

| | `cerny2007_extraction.md` (Cerny & Briffod) | **this paper** | `zhao2019_extraction.md` |
|---|---|---|---|
| system | [¹³C₅]-xylose + Cys + thiamine, pH 4-7, 145 °C, 20 min | 1 : 1 [¹³C₅]/[¹²C]-xylose + Cys + thiamine ± peptide, pH 4.8, 140 °C, 90 min | [¹³C₅]-TTCA + Cys + Gly + Gly-Amadori, pH 5.5, 140 °C, 30 min |
| FFT carbon from sugar | **> 99 %** | **none detected** | **77 %** (system b) |
| MFT carbon from sugar | 10-25 % | none detected | 69 % (system b) |
| quantity | peak areas | "relative concentration" | relative peak areas |

## 2. Methods as they matter to a model

- **Charge (Table 1, mg in 25 mL distilled water).** Xyl 50 + [¹³C₅]-Xyl 50 (98 % enrichment) +
  thiamine ("VB1") 10 + Cys 57 + one peptide 57 (A-D); peptide + xylose only (E-G). → total xylose
  **26.2 mmol/L** (13.3 unlabelled + 12.9 labelled), Cys **18.8 mmol/L**, thiamine **1.19 mmol/L**
  if the hydrochloride (1.51 if free base; form not stated), peptides **~2.4 / 2.3 / 2.0 mmol/L**
  (peptide-1 943 Da, peptide-2 1007 Da, peptide-3 1144 Da by average residue masses — mine).
- **Peptides** (Apeptides, Shanghai, "95% enrichment" [sic, purity]): peptide-1
  **Asp-Prp-Glu-Phe-Pro-Asp-Pro-Lys** [sic; "Prp" presumably Pro]; peptide-2
  **Gly-Val-Thr-Val-Phe-Glu-Asp-Leu-Lys**; peptide-3 **Glu-Thr-Asp-Asn-Leu-Asp-Ser-Val-Pro-Arg**.
- **Reaction.** "adjusted to **pH 4.8** and heated at **140 °C for 90 min** in a high-pressure stainless
  reactor (Parr Instrument Company) with stirred at **350 rpm**"; cooled with running tap water. Acid
  or base used for the adjustment not stated; reactor volume (hence headspace) not stated;
  atmosphere not stated; **final pH not measured**.
- **Volatiles.** Purge and trap (Stratum, Tekmar): 10 mL reaction mixture + 1 µL internal standard
  **2-methyl-3-heptanone 1.68 µg/µL** (in methanol per 2.3; "in ethyl ether" per 2.6), 50 °C 30 min,
  N₂ 40 mL/min purge 12 min at 50 °C, desorb 250 °C 2 min, bake 280 °C 5 min. GC-MS 7890A / Triple
  Quad 7000B, **DB-WAX UI** 30 m × 0.25 mm × 0.25 µm, He 1.2 mL/min, 40 °C 3 min → 5 °C/min → 200 °C →
  15 °C/min → 260 °C 5 min; Sniffer 9000 olfactometer, three panellists. Identification by NIST 2.0,
  LRI (C7-C22) and odour.
- **Quantification.** "concentration of an odorant = (concentration of internal standard × peak area
  of compound × response factor) / peak area of the internal standard; since the peak area ratio was
  equal to concentration ratio" — i.e. **response factor = 1**; reported as "Relative concentration"
  with **no unit**; "mean value of triplicates", **no SD printed**.
- **CAMOLA.** 1 : 1 fully labelled : unlabelled D-xylose; isotopomer ratios from relative signal
  intensities [M⁺] to [M⁺+5] versus the unlabelled trial.
- **Non-volatiles.** Agilent 6530 ESI-Q-TOF + 1290 HPLC, Zorbax Eclipse XDB-C18 2.1 × 50 mm 1.8 µm,
  25 °C, 5 → 30 % MeOH/H₂O (0.1 % TFA) over ~30 min, 0.2 mL/min; positive ESI, m/z 100-3000.

## 3. Tables re-typed

### Table 1 (verbatim). Model reactions — Amount (mg), each in 25 mL water

| | A | B | C | D | E | F | G |
|---|---|---|---|---|---|---|---|
| Xyl | 50 | 50 | 50 | 50 | 50 | 50 | 50 |
| [¹³C₅]Xyl | 50 | 50 | 50 | 50 | 50 | 50 | 50 |
| Thiamine | 10 | 10 | 10 | 10 | | | |
| Cys | 57 | 57 | 57 | 57 | | | |
| Peptide-1 | | 57 | | | 57 | | |
| Peptide-2 | | | 57 | | | 57 | |
| Peptide-3 | | | | 57 | | | 57 |

### Table 2 (verbatim). Volatile compounds identified in model reaction — "Relative concentration" (unitless as printed)

"—" as printed; an empty cell is printed as empty in the source.

| Compound | LRI | Odor | A | B | C | D | E | F | G |
|---|---|---|---|---|---|---|---|---|---|
| *Thiophene* | | | | | | | | | |
| Thiophene | 1029 | Pungent | 8.384 | 4.215 | 4.429 | 8.654 | — | — | — |
| 2-Methyl-thiophene | 1094 | Onion | 6.991 | 6.398 | 2.659 | 4.419 | — | — | — |
| 2-Ethyl-thiophene | 1165 | Styrene | — | 0.122 | — | 0.217 | — | — | — |
| 2,4-Dimethyl-thiophene | 1205 | Metallic, sulfurous | — | 0.223 | 0.352 | 0.236 | — | — | — |
| 3-Thiophenecarboxaldehyde | 1687 | Sulfurous | 0.683 | 0.357 | 0.687 | 0.216 | — | — | — |
| 2-Thiophenecarboxaldehyde | 1702 | Sulfurous | 0.846 | 0.530 | 0.402 | 0.455 | — | — | — |
| 5-Methyl-2-thiophenecarboxaldehyde | 1721 | Sulfurous | 0.778 | 0.793 | 0.815 | 1.396 | — | — | — |
| 2-Methoxy-5-methyl-thiophene | 1776 | Fatty | — | 1.110 | 1.162 | 2.349 | — | — | — |
| Thieno[3,2-b]thiophene | 1887 | Metallic | 1.259 | 0.386 | 0.344 | 0.476 | — | — | — |
| 2,5-Thiophenedicarboxaldehyde | 1924 | Sulfurous | — | — | — | 0.250 | — | — | — |
| 2-Methyl-thieno[2,3-b]thiophene | 1955 | | 0.256 | 0.088 | 0.064 | 0.081 | — | — | — |
| 3-(Vinylthio)thiophene | 2029 | | — | 0.282 | 0.092 | 0.199 | — | — | — |
| *Thiazole* | | | | | | | | | |
| Thiazole | 1257 | Meaty | 0.837 | 0.672 | 0.540 | 0.645 | — | — | — |
| 4-Methyl-5-thiazoleethanol | 2309 | Roasted | — | — | 0.199 | — | — | — | — |
| *Thiopenone* [sic] | | | | | | | | | |
| Dihydro-2-methyl-3(2H)-thiophenone | 1535 | Meaty | 0.818 | 2.251 | 1.572 | 2.313 | — | — | — |
| *Furans and their derivatives* | | | | | | | | | |
| 2,3-Dihydrofuran | 1049 | | 2.443 | 1.826 | 1.582 | 3.452 | — | — | — |
| 2-Butyl-furan | 1065 | Noncharacteristic | 0.759 | — | — | 0.446 | — | — | — |
| 2-Pentyl-furan | 1213 | Fatty | — | — | — | — | — | 0.240 | — |
| Dihydro-2-methyl-3(2H)-Furanone | 1273 | Fruity | 2.112 | 0.406 | 0.903 | 0.825 | — | — | — |
| **2-Methyl-3-furanthiol** | 1334 | Meaty | **0.042** | **0.054** | **0.043** | **0.243** | — | — | — |
| **2-Furfurylthiol** | 1436 | Roasted, meaty | **1.133** | **1.717** | **0.294** | **1.015** | — | — | — |
| **Furfural** | 1472 | Caramel | **142.4** | **148.2** | **168.3** | **152.4** | **113.7** | **59.73** | **160.3** |
| 1-(2-Furanylmethyl)-1H-pyrrole | 1834 | Green, hay | — | 0.067 | 0.090 | 0.121 | 0.184 | 0.078 | 0.117 |
| *Aldehydes, ketones* | | | | | | | | | |
| 3-Methyl-butanal | 897 | malty Peach, fatty [sic] | — | — | — | 15.411 | 19.700 | 11.825 | 28.606 |
| 2-Pentanone | 1012 | Fruity, sweet | 2.711 | 1.069 | 3.399 | 2.733 | 10.710 | 9.299 | 8.202 |
| 2,3-Pentanedione | 1067 | Fruity | 2.705 | 2.337 | 1.664 | 1.713 | — | — | — |
| Hexanal | 1083 | Green, grass | — | — | 0.571 | 0.591 | 13.262 | 8.463 | 7.986 |
| 5-Methyl-2-hexanone | 1178 | Sweet | — | — | — | — | 11.267 | 6.128 | 7.782 |
| Heptanal | 1182 | Green | — | — | — | 0.426 | — | — | — |
| 4-Octanone | 1217 | Rusty | — | — | — | — | — | 0.506 | 0.623 |
| Octanal | 1280 | Fruity, sour | 0.652 | — | — | 0.205 | | 2.065 | 3.949 |
| 2-Nonanone | 1373 | Musty | — | — | — | — | | 0.274 | — |
| Nonanal | 1387 | Green, stale | 1.002 | 0.105 | — | 0.461 | | 1.897 | 3.318 |
| (E)-2-Octenal | 1426 | Green nutty | — | — | — | — | 1.059 | 0.608 | 0.730 |
| Undecanal | 1433 | Fatty | — | — | — | — | | 0.326 | — |
| Benzaldehyde | 1522 | Almond | 0.968 | 3.439 | 2.995 | 0.648 | 8.897 | 7.064 | 2.428 |
| (E)-2-Nonenal | 1532 | Cucumber | — | — | — | — | 0.428 | 0.390 | 0.625 |
| 4-Methyl-benzaldehyde | 1649 | | — | — | 0.259 | 0.277 | 0.582 | 0.519 | 0.772 |
| 3-Ethyl-benzaldehyde | 1706 | | — | — | — | — | — | 0.090 | 0.128 |
| *Esters* | | | | | | | | | |
| Heptyl propionate | 1036 | | — | — | — | — | 1.862 | — | 0.901 |
| Ethyl hexanoate | 1221 | Sweet | — | — | — | — | — | 0.540 | 1.059 |
| Methyl-2-ethyl hexanoate | 1247 | | — | — | — | — | 1.133 | 0.387 | 0.381 |
| *Miscellaneous* | | | | | | | | | |
| 3-Ethyl-cyclopentanone | 1331 | | | — | — | — | — | 0.729 | — |
| 2-Butyl-1-octanol | 1375 | | | — | — | — | — | — | 0.849 |
| 2-Methyl-2-cyclopenten-1-one | 1377 | | | — | — | 0.078 | — | | — |
| 3,5-Octadien-2-ol | 1404 | | | — | — | — | — | 0.267 | — |
| 5-Ethylcyclopent-1-enecarboxaldehyde | 1413 | | | — | — | — | 1.353 | 0.445 | 0.528 |
| Heptanol | 1449 | | | — | — | — | 2.843 | 1.984 | 1.816 |
| Camphor | 1515 | | | — | — | — | 3.087 | 2.219 | 2.428 |
| 3,7-Dimethyl-1,6-octadien-3-ol | 1538 | | | — | — | — | 0.563 | — | |
| 3-Butyl-cyclopentanone | 1541 | | | — | — | — | — | 0.328 | 0.446 |
| Octanol | 1550 | | | — | — | — | 2.308 | 1.858 | 2.452 |

Footnotes (verbatim in substance): ᵃ identified by mass spectra and retention indices; ᵇ RIs vs
n-alkanes C7-C22 as external references; ᶜ odour descriptions at the Sniffer 9000. The odour column
is offset by one row in places in the source (e.g. "Styrene" for 2-ethylthiophene); reproduced as
printed.

### Table 3 (verbatim). Proportion of isotopomers in compounds formed from the model reactions (%)

(a) systems A-D (pooled, not stated which):

| No. | Compound | m/z | [M⁺] | [M⁺+1] | [M⁺+2] | [M⁺+3] | [M⁺+4] | [M⁺+5] |
|---|---|---|---|---|---|---|---|---|
| 1 | 2,3-Dihydrofuran | 70 | 54 | | | | 46 | |
| 2 | Thiophene | 84 | 77 | | | | 23 | |
| 3 | Furfural | 96 | 56 | | | | | 44 |
| 4 | 2-Butyl-furan | 124 | 56 | | | | | 44 |
| 5 | 2-Thiophenecarboxaldehyde | 112 | 53 | 7 | 8 | 5 | 14 | 13 |
| 6 | 5-Methyl-2-thiophenecarboxaldehyde | 126 | 51 | | | 26 | 23 | |
| 7 | 2-Methoxy-5-methyl-thiophene | 128 | 60 | | | | 40 | |

(b) Proportion of isotopomers formed from the thermal reaction of E, F and G:

| No. | Compound | m/z | [M⁺] | [M⁺+2] | [M⁺+3] | [M⁺+4] | [M⁺+5] | [M⁺+9] |
|---|---|---|---|---|---|---|---|---|
| 1 | Furfural | 96 | 56 | | | | 44 | |
| 2 | 1-(2-Furanylmethyl)-1H-pyrrole | 147 | 43 | | | | 42 | 15 |

**MFT and FFT do not appear in Table 3.** The text says they "did not appear labeled molecules".

### Arithmetic (mine)

**1. Expected label share.** With 13.3 mM unlabelled and 12.9 mM labelled xylose, a purely
xylose-derived C5 product should be ~49 % M+5 (ignoring isotope effects). Furfural 44 % M+5 and
2-butylfuran 44 % M+5 are close (≈ 0.9 of expected). Thiophene 23 % M+4 means **~47 % of thiophene
is xylose-derived** and the remainder from an unlabelled source (cysteine/peptide/thiamine).

**2. MFT is a trace product here.** 0.042-0.243 relative units against 142-168 for furfural —
3-4 orders of magnitude lower, so an M+5 cluster at half of an already minor ion would sit at the
noise floor. The "no label" claim for MFT and FFT is therefore weak evidence.

**3. Peptides lower some sulfur volatiles relative to A** (thiophene, 2-methylthiophene, thieno-
thiophene) and raise others (2-methoxy-5-methylthiophene appears only with peptides). FFT: B / A =
1.52, C / A = 0.26, D / A = 0.90. MFT: D / A = 5.8, others ≈ 1. Single replicate means without SD —
treat as direction only.

## 4. Kinetic numbers the repository can use

**None.** No rate, no time series, no intermediate concentration, no unit on the volatile amounts.

| step / observable | quantity | value | unit | conditions | order | source | evidence class |
|---|---|---|---|---|---|---|---|
| furfural is xylose-derived (intact C5) | M⁺ / M⁺+5 | 56 / 44 | % | 1 : 1 labelled xylose, 140 °C 90 min, pH 4.8 | — | Table 3a, 3b | within_study_ratio |
| thiophene from C4 of xylose | M⁺ / M⁺+4 | 77 / 23 | % | same | — | Table 3a | within_study_ratio — **out of the module's scope** (thiophenes) |
| MFT, FFT carbon origin | "no labelled molecules" | — | — | same | — | text only | **threshold, contradicted by corpus** — do not use |
| volatile amounts | Table 2 | unitless | — | — | — | Table 2 | **level_only** (no unit, RF = 1, no SD) |
| pyrazines | not detected | — | — | pH 4.8, 140 °C | — | text | threshold (negative) |

## 5. What this can and cannot anchor in the model

**Can:** nothing quantitative. At most: (i) furfural from intact xylose under Cys + thiamine at
pH 4.8 (consistent with `r_tdp_fur`), (ii) a qualitative "no pyrazines at pH 4.8" negative.

**Cannot:** anything on TTCA, Amadori compounds, cysteine, H2S, MFT/FFT rates or yields. The pots
contain thiamine (A-D) or no sulfur at all (E-G), and peptides with no engine species. The MFT/FFT
carbon-origin claim conflicts with Cerny & Briffod 2007 and Zhao 2019 and with this paper's own
furfural labelling.

## 6. Flags

1. **The FFT "no label" claim contradicts the corpus and the paper itself.** FFT forms from furfural
   + H2S; this paper's furfural is 44 % M+5; Cerny & Briffod 2007 measure FFT > 99 % ¹³C₅ from
   [¹³C₅]-xylose in the same ternary system; Zhao 2019 measure 77 % from [¹³C₅]-TTCA. Probably a
   detection-limit artefact on a trace signal (arithmetic 2).
2. **Abstract and text disagree on the peptide cleavage site.** Abstract: peptide-2 "capable of
   splitting in the position of aspartyl and **prolyl**"; text: "in the position of aspartyl and
   **leucyl**". Peptide-2 contains no proline.
3. **Typos in the materials:** peptide-1 "Asp-**Prp**-Glu-...", "95% enrichment" for peptide purity,
   "Thiopenone" as a class heading.
4. **Internal standard solvent stated two ways** (methanol in sec. 2.3, ethyl ether in sec. 2.6).
5. **"Relative concentration" has no unit**, response factors are set to 1, and no SD is printed
   despite "mean value of triplicates".
6. **Table 3a does not say which of A-D it pools**, and omits MFT/FFT entirely.
7. **pH set once (4.8, agent not stated), never re-measured**; reactor volume / headspace unknown;
   atmosphere not stated.
8. **Thiamine form (free base vs hydrochloride) not stated** — 1.19 vs 1.51 mmol/L.
9. **Different BTBU group** (Song, Laboratory of Molecular Sensory Science) from Cao 2017 / Zhao
   2019 (Xie group); the protocol parent Lv & Song 2010 is not on disk.
10. **Registry:** `2_methyl_3_furanthiol`, `2_furfurylthiol`, `2_methylthiophene`, `furfural`,
    `3_methylbutanal`, `hexanal`, `benzaldehyde`, `2_pentylfuran` present; thiamine, cysteine,
    xylose and the peptides absent. Paper not in `papers.yml`.

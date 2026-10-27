# Cerny 2007b — EXTRACTION
### *"Origin of carbons in sulfur-containing aroma compounds from the Maillard reaction of xylose, cysteine and thiamine"*
### [13C5]xylose + cysteine + thiamine at 300 : 100 : 100 mM, pH 5.00 (0.5 M K-phosphate), 145 °C / 20 min, plus the two single-route controls. It gives carbon origin for 10 sulfur volatiles. MFT is 46 % xylose-derived, 3-mercapto-2-pentanone 40 %, FFT and furfural >99 %. With no cysteine, everything except FFT and furfural comes from thiamine; with no thiamine, everything is xylose. Isotope ratios only: no amounts of any kind. Read 2026-10-27.

**Full read of `data/articles/cerny2007b.pdf`** (7 pages, the Elsevier "ARTICLE IN PRESS" PDF that
carries the final volume and pages; clean text layer). Tables 1–4 were checked against a 130-dpi
raster of pp. 1310–1312; the text layer prints "o" for "<" and "4" for ">", corrected below.
Companion: `cerny2007_extraction.md` (Cerny & Briffod 2007 JAFC, the same chemistry at half the
concentration across pH 4–7, which cites this paper as ref 9). **Note:** `k3_final_parameter_inventory.md`
(l. 1105) quotes a `cerny2007b_extraction.md` that was never on disk or in git. This file is the
first actual primary-source read; the k3 numbers it quoted are confirmed in §6.

**Provenance:** **[M]** measured, printed as a number · **[Q]** qualitative · **[Z]** derived by me
from the paper's own numbers, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Author | Christoph Cerny |
| Title | "Origin of carbons in sulfur-containing aroma compounds from the Maillard reaction of xylose, cysteine and thiamine" |
| Venue | *LWT* **40** (2007) **1309–1315** |
| DOI | 10.1016/j.lwt.2006.09.008 |
| Affiliation | Firmenich SA, Corporate Research, Meyrin 2, Geneva, Switzerland |
| Dates | received 17 Jul 2006, revised 19 Sep 2006, accepted 20 Sep 2006 |
| PDF character | Table 1 (recipes), Table 2 (full molecular-ion clusters for 4 thiols), Tables 3/5/6 (isotopomer splits per system), Table 4 (one thiamine product); Fig. 1 structures, Fig. 2 scheme. **No peak areas, no concentrations, no replicates stated** |

---

## 1. THE ONE-PARAGRAPH ANSWER

Cerny heats xylose 300 mM, cysteine 100 mM and thiamine HCl 100 mM in 1.00 mL of 0.5 M K-phosphate
at pH 5.00, 145 °C, 20 min, in septum-closed 2 mL Teflon vials. Each pot is run once with
[12C] and once with [13C5]xylose (98 %), then the cysteine-free and thiamine-free controls are
run the same way. He reads the M / M+5 isotopomer ratio of each sulfur volatile by HS-SPME-GC-MS.
**In the full ternary, MFT is 54 : 46 and 3-mercapto-2-pentanone 60 : 40 (unlabelled : [13C5]),
so ~half of each skeleton comes from xylose. FFT and furfural are >99 % xylose.**
4,5-Dihydro-2-methyl-3(2H)-furanone, 3-mercapto-2-butanone, 4,5-dihydro-2-methyl-3-furanthiol and
the thiophenone are ≥95 % thiamine. 2-Methylthiophene is 62 : 38. 2-Methyl-3-(methylthio)furan is
a four-isotopomer mix in which **70 % of the S-methyl carbons are xylose-derived** [Z].
**Without cysteine**, MFT, 3-mercapto-2-pentanone, 2-methylthiophene and 2,3-pentanedione are
>99 % thiamine. **FFT is still >99 % xylose**, so thiamine supplied the sulfur to a xylose skeleton
at pH 5. **Without thiamine**, MFT and 3-mercapto-2-pentanone are >95 % xylose, and
2-methyl-3-(methylthio)furan still forms, with its methyl 64 % from xylose. Cysteine therefore
acts as the sulfur donor that unlocks the sugar route, not as a carbon source.

---

## 2. METHODS [M] (§2, p. 1310)

| item | as printed | [Z] / note |
|---|---|---|
| Charges | cysteine **12.1 mg**, xylose **45.0 mg** (or [13C5]xylose 45.0 mg), thiamine hydrochloride **33.7 mg**, in **1.00 mL** | Cys **99.9 mM**; xylose **299.7 mM** ([13C5]: 45.0/155.09 = **290.2 mM**, 3 % fewer moles); thiamine **99.9 mM** |
| Molar ratio | "cysteine, xylose and thiamine was 1:3:1" | |
| Buffer | **K-phosphate, 0.5 mol/L, pH 5.00** (§2.2 and Table 1 fn a) | pH set by the buffer at room temperature; **post-reaction pH not reported** |
| Vessel | **2 mL septum-closed Teflon vial** (Infochroma), stirred, heated metal block (Reacti-Therm, Pierce) | liquid 1.00 mL, so **~1 mL headspace**, atmosphere not stated (presumably air) |
| Temperature / time | **145 °C, 20 min** | come-up time not stated |
| Leak check | "Vials were weighed before and after the reaction, and no leakages were detected" | |
| Controls | "Similar single experiments were carried out omitting cysteine and thiamine, respectively" | **[!] single experiments, unreplicated**; the ternary replication is not stated either |
| Analysis | HS-SPME, **50/30 µm DVB/CAR on PDMS**, **35 min at 40 °C** over the sample in **10 mL glass vials**; desorption 3.5 min at 250 °C (0.75 mm liner); Agilent 6890A/5973; HP-5MS 30 m × 0.25 mm × 0.25 µm; 40 °C (5 min) → 5 °C/min → 260 °C → 15 °C/min → 280 °C; EI 70 eV, m/z 15–350 | the delay between reaction and analysis is not stated |
| Quantity reported | **isotope ratios from molecular-ion (or M−1) intensities only** | SPME partition is the same for isotopomers, so ratios are robust to the fibre; **absolute amounts are not reported at all** |
| Identification | authentic references for all except 8 and 12 (tentative, literature MS/RI) | |

---

## 3. TABLES AS PRINTED

### Table 1. "Model reactions" (fn a: "Reaction in phosphate buffer (0.5 mol/l; pH 5.00; 1.00 ml) at 145 °C (20 min)."). Amount (mg)

| | A | B | C | D | E | F |
|---|---:|---:|---:|---:|---:|---:|
| Cysteine | 12.1 | 12.1 | | | 12.1 | 12.1 |
| Xylose | 45.0 | | 45.0 | | 45.0 | |
| [13C5]xylose | | 45.0 | | 45.0 | | 45.0 |
| Thiamine HCl | 33.7 | 33.7 | 33.7 | 33.7 | | |

So: A/B = ternary (unlabelled/labelled), C/D = **no cysteine** (xylose + thiamine), E/F = **no
thiamine** (xylose + cysteine).

### Table 2. "Isotope ratios (%) of selected sulfur volatiles from the reaction of thiamine and cysteine with xylose (Xyl) or [13C5]xylose ([13C5]Xyl), respectively" (ternary A vs B)

2-Methyl-3-furanthiol (5) and 3-mercapto-2-pentanone (6):

| m/z | MFT Xyl | MFT [13C5]Xyl | | m/z | MP Xyl | MP [13C5]Xyl |
|---:|---:|---:|---|---:|---:|---:|
| 112 | 0.7 | 0.3 | | 116 | 8.8 | 4.6 |
| 113 | 17.7 | 9.6 | | 117 | <0.1 | <0.1 |
| **114** | **72.4** | **38.9** | | **118** | **81.7** | **48.7** |
| 115 | 5.4 | 2.9 | | 119 | 5.2 | 3.2 |
| 116 | 3.6 | 2.0 | | 120 | 4.1 | 2.6 |
| 117 | 0.2 | 1.0 | | 121 | 0.2 | 3.6 |
| 118 | <0.1 | 10.3 | | 122 | <0.1 | 2.3 |
| **119** | <0.1 | **32.6** | | **123** | <0.1 | **33.0** |
| 120 | <0.1 | 0.8 | | 124 | <0.1 | 0.4 |
| 121 | <0.1 | 1.6 | | 125 | <0.1 | 1.6 |

2-Furfurylthiol (7) and 4,5-dihydro-2-methyl-3-furanthiol (8):

| m/z | FFT Xyl | FFT [13C5]Xyl | | m/z | DHMFT Xyl | DHMFT [13C5]Xyl |
|---:|---:|---:|---|---:|---:|---:|
| 112 | 0.7 | <0.1 | | 114 | 4.7 | 2.9 |
| 113 | <0.1 | <0.1 | | 115 | 7.4 | 6.9 |
| **114** | **89.3** | 0.3 | | **116** | **78.5** | **80.4** |
| 115 | 5.4 | 0.1 | | 117 | 5.1 | 5.3 |
| 116 | 4.4 | 1.5 | | 118 | 4.0 | 4.0 |
| 117 | 0.2 | 0.9 | | 119 | 0.3 | 0.3 |
| 118 | <0.1 | 7.4 | | 120 | <0.1 | <0.1 |
| **119** | <0.1 | **84.7** | | 121 | <0.1 | <0.1 |
| 120 | <0.1 | 1.1 | | 122 | <0.1 | <0.1 |
| 121 | <0.1 | 4.0 | | 123 | <0.1 | <0.1 |

### Table 3. "Proportion of isotopomers formed from [13C5]xylose, cysteine and thiamine" (ternary, B)

Footnotes: a, identified vs authentic references unless marked e; b, RI on HP-5MS; c, "Values represent
the proportion between labeled and the indicated unlabeled isotopomers and are based on the abundance
of the respective analysed ion signals"; d, "Calculations were based on the (M+ − 1) ion"; e, tentative.

| No. | compound | RI | m/z (M+) | m/z analysed | unlabeled C (%) | labeled C (%) | no. of labeled C |
|---:|---|---:|---:|---|---:|---:|---|
| 1 | 2-methylthiophene (d) | 769 | 97 (d) | 97; 102 | 62 | 38 | 5 |
| 2 | 4,5-dihydro-2-methyl-3(2H)-furanone | 800 | 100 | 100; 105 | >99 | <1 | 5 |
| 3 | 3-mercapto-2-butanone | 807 | 104 | 104; 108 | >95 | <5 | 4 |
| 4 | 2-furaldehyde | 831 | 96 | 96; 101 | <1 | >99 | 5 |
| 5 | **2-methyl-3-furanthiol** | 871 | 114 | 114; 119 | **54** | **46** | 5 |
| 6 | **3-mercapto-2-pentanone** | 901 | 118 | 118; 123 | **60** | **40** | 5 |
| 7 | **2-furfurylthiol** | 915 | 114 | 114; 119 | **<1** | **>99** | 5 |
| 8 | 4,5-dihydro-2-methyl-3-furanthiol (e) | 939 | 116 | 116; 121 | >99 | <1 | 5 |
| 9 | 2-methyl-3-(methylthio)furan | 950 | 128 | 128; 129; 133; 134 | 6 | 22; 24; 48 | 1; 5; 6 |
| 10 | 4,5-dihydro-2-methyl-3(2H)-thiophenone | 988 | 116 | 116; 121 | >95 | <5 | 5 |

### Table 4. "Isotope ratios (%) of 1-methylbicyclo[3.3.0]-4-thia-2,8-dioxooctane from the reactions thiamine/xylose (TX), thiamine/[13C5]xylose (T[13C5]X), thiamine/cysteine/xylose (TCX) and thiamine/cysteine/[13C5]xylose (TC[13C5]X)"

| m/z | TX | T[13C5]X | TCX | TC[13C5]X |
|---:|---:|---:|---:|---:|
| 145 | <0.1 | <0.1 | <0.1 | <0.1 |
| 146 | **88.9** | **87.0** | **86.1** | **72.4** |
| 147 | 6.3 | 8.0 | 6.1 | **17.9** |
| 148 | 4.5 | 4.5 | 7.0 | 4.7 |
| 149 | 0.3 | 0.3 | 0.5 | 0.9 |
| 150 | <0.1 | <0.1 | 0.2 | 0.2 |
| 151 | <0.1 | <0.1 | <0.1 | 0.1 |
| 152 | <0.1 | <0.1 | <0.1 | <0.1 |

Text: "Apparently one xylose carbon is integrated to about 15% into the molecule" (via [13C]formaldehyde,
presumably C-3). [Z] (17.9 − 6.1) / 90.3 ≈ 13 %.

### Table 5. "Proportion of isotopomers formed from [13C5]xylose and thiamine" (NO CYSTEINE, D)

| No. | compound | RI | m/z (M+) | m/z analysed | unlabeled C (%) | labeled C (%) | no. of labeled C |
|---:|---|---:|---:|---|---:|---:|---:|
| 11 | 2,3-pentanedione | 696 | 100 | 100; 105 | >99 | <1 | 5 |
| 1 | 2-methylthiophene (d) | 769 | 97 (d) | 97; 102 | >99 | <1 | 5 |
| 2 | 4,5-dihydro-2-methyl-3(2H)-furanone | 800 | 100 | 100; 105 | >99 | <1 | 5 |
| 4 | 2-furaldehyde | 831 | 96 | 96; 101 | <1 | >99 | 5 |
| 5 | **2-methyl-3-furanthiol** | 871 | 114 | 114; 119 | **>99** | **<1** | 5 |
| 6 | **3-mercapto-2-pentanone** | 901 | 118 | 118; 123 | **>99** | **<1** | 5 |
| 7 | **2-furfurylthiol** | 915 | 114 | 114; 119 | **<1** | **>99** | 5 |
| 8 | 4,5-dihydro-2-methyl-3-furanthiol (e) | 939 | 116 | 116; 121 | >99 | <1 | 5 |

Text (§3.2): the two main peaks of the xylose + thiamine chromatogram are compound 12 and furfural (4).

### Table 6. "Proportion of isotopomers formed from [13C5]xylose and cysteine" (NO THIAMINE, F)

| No. | compound | RI | m/z | m/z analysed | unlabeled C (%) | labeled C (%) | no. of labeled C |
|---:|---|---:|---:|---|---:|---:|---|
| 1 | 2-methylthiophene (d) | 769 | 97 | 97; 102 | 46 | 54 | 5 |
| 4 | 2-furaldehyde | 831 | 96 | 96; 101 | <1 | >99 | 5 |
| 5 | **2-methyl-3-furanthiol** | 871 | 114 | 114; 119 | **<5** | **>95** | 5 |
| 6 | **3-mercapto-2-pentanone** | 901 | 118 | 118; 123 | **<5** | **>95** | 5 |
| 7 | 2-furfurylthiol | 915 | 114 | 114; 119 | <1 | >99 | 5 |
| 9 | 2-methyl-3-(methylthio)furan | 950 | 128 | 128; 129; 133; 134 | <1 | 36; 64 | 5; 6 |

Text (§3.3): compounds 2 and 8 (thiamine markers) not detected; 3-mercapto-2-butanone (3) not detected.

---

## 4. CARBON-ORIGIN MATRIX (% of the skeleton from xylose) [M], with [Z] cross-checks

| compound | ternary (100/100/300 mM) | no cysteine | no thiamine | origin verdict |
|---|---:|---:|---:|---|
| furfural | >99 | >99 | >99 | xylose, intact C5 |
| **2-furfurylthiol** | **>99** | **>99** | **>99** | xylose skeleton; S from cysteine **or thiamine** |
| **2-methyl-3-furanthiol** | **46** | **<1** | **>95** | two routes: thiamine (C + S) and xylose (C) + cysteine (S) |
| **3-mercapto-2-pentanone** | **40** | **<1** | **>95** | same two routes |
| 2-methylthiophene | 38 | <1 | 54 (46 % unlabelled) | thiamine, xylose, and a non-xylose route with no thiamine (cysteine carbons, per the author) [!] |
| 2-methyl-3-(methylthio)furan | furan 72, S-CH3 70 [Z] | — | furan >99, S-CH3 64 | furan follows the xylose share (higher than MFT's own 46 %); the methyl is mostly xylose C1 |
| 4,5-dihydro-2-methyl-3(2H)-furanone | <1 | <1 | n.d. | thiamine |
| 3-mercapto-2-butanone | <5 | — | n.d. | thiamine |
| 4,5-dihydro-2-methyl-3-furanthiol | <1 | <1 | n.d. | thiamine (via 5-hydroxy-3-mercapto-2-pentanone, Fig. 2) |
| 4,5-dihydro-2-methyl-3(2H)-thiophenone | <5 | — | — | thiamine |
| 2,3-pentanedione | — | <1 | — | thiamine (not xylose) at pH 5 with no cysteine |
| bicyclic 12 | ~13–15 % singly labelled | ~2 % [Z] | — | thiamine, with one xylose C (formaldehyde) when cysteine is present |

**[Z] Cross-checks against Table 2:** MFT [13C5] cluster m/z 118–121 = 45.3 % vs 112–117 = 54.7 %,
reproducing 54 : 46. MP 121–125 = 40.9 % vs 116–120 = 59.1 %, reproducing 60 : 40. FFT ≥ 99 % at
m/z 116–121. DHMFT patterns are identical with and without label. The table is internally
consistent.

**[Z] Methylthio-MFT isotopomers (ternary):** 128 (6 %) = both parts unlabelled; 129 (22 %) = S-13CH3
on an unlabelled furan; 133 (24 %) = [13C5]furan + 12CH3; 134 (48 %) = both labelled. Furan from
xylose 24 + 48 = **72 %**; methyl from xylose 22 + 48 = **70 %**.

**[!] Prose errors (already logged in k3 l. 1609):** the Discussion says "both precursors had a
similar share in the carbon skeleton of MFT (54% 5-times labeled) and MP (60% 5-times labeled, cf.
Table 5)". The tables give 46 % and 40 % labelled, and the ternary numbers are in Table 3, not
Table 5. **Ingest the tables, not the prose.** The abstract's "4,5-dihydro-2-methyl-3(2H)-thiophene"
means the thiophenone (10).

**Literature numbers quoted in the paper ([cited], not measured here):** Bolton et al. 1994,
[34S]cysteine + thiamine + xylose: only 8 % of the MFT carried 34S (thiamine was the dominant S
donor), and no MFT at all without thiamine. Hincelin et al. 1992: MFT from heated thiamine rose
4–5× with xylose (pH 7, 100 °C, 1 h, 20-fold xylose). Grosch / Zeiler-Hilgart: thiamine far more
efficient than cysteine/ribose at >100× lower concentrations.

---

## 5. RELATION TO CERNY & BRIFFOD 2007 (`cerny2007_extraction.md`)

Same lab, same chemistry (0.5 M K-phosphate, 145 °C / 20 min, Teflon, HS-SPME DVB/CAR/PDMS, same
GC). C&B run at **half the concentration** (150 / 50 / 50 mM) across pH 4.0–7.0 with peak areas.
This paper is **"A" in C&B's Table 5**.

| quantity at pH 5 | this paper (300/100/100 mM) | C&B (150/50/50 mM) | ratio |
|---|---:|---:|---:|
| MFT xylose share | **46 %** | 15 % | 3.1× |
| 3-mercapto-2-pentanone xylose share | **40 %** | 14 % | 2.9× |
| FFT xylose share | >99 % | >99 % | — |
| furfural | >99 % | >99 % | — |
| 4,5-dihydro-2-methyl-3-furanthiol | <1 % | <1 % | — |
| methylthio-MFT isotopomers 128 / 129 / 133 / 134 | 6 / 22 / 24 / 48 | 39 / 3 / 8 / 50 | the furan's xylose share 72 → 58 %; the methyl's 70 → 53 % |
| 2-mercapto-3-pentanone | not reported | absent at pH 4–5.5, 94–95 % xylose at pH 6–7 | — |

**Adding Wang, Yang & Song 2012** (`wang2012_extraction.md`): thiamine 1.2 mM, xylose 26 mM,
Cys ~10 mM, unbuffered pH ~4.8, 132 °C / 90 min, 1:1 label mix. **MFT <~2 % xylose.** Under
different T, t, ratio and buffer, this extends the sign of the series: **the xylose share of MFT
falls monotonically as the precursors are diluted (46 → 15 → <2 %)**. Cerny's own explanation
is that "in diluted systems thiamine is relatively more efficient". The model's structure gives
the same direction: thiamine → HMP is first order, the sugar route is second order in
(osone × sulfide).

---

## 6. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Already in the repo** (FIT_HOLDOUT_DECLARATION l. 140–141; k3 §(ii); `generate_kinetic_core_b2_holdout.py`
§11): the ternary 54 : 46 = **FIT**; the two controls (no-Cys MFT ≥ 99 % thiamine; no-THI MFT
≥ 95 % sugar) = **HOLD-OUT**. This read **confirms all three numbers** against the tables. Two
small corrections to how they are carried:
- the labelled pots B/D/F contain **290 mM** [13C5]xylose, not 300 (the generator uses
  `PENT: 299.7`). The 3 % difference is negligible but should be stated;
- k3's "system A" is C&B Table 5's label for this paper's conditions. In this paper's own
  Table 1, A is the **unlabelled** ternary and the split comes from **B**. The conditions are
  identical, so only the naming differs.

**New things it can anchor or test (none of them amounts):**

1. **The 3-mercapto-2-pentanone split, a second FIT-grade branching ratio alongside MFT:** 60 : 40
   thiamine : xylose in the ternary (Table 3), with controls >99 : <1 (no Cys) and <5 : >95 (no THI).
   **[!] STRUCTURAL CONTRADICTION.** The model's `MP2P` (3-mercapto-2-pentanone, labelled the
   "THIAMINE-diagnostic isomer" in `species_sulfur.py`) has **one producer, `r_hmp_mp2p`, from
   thiamine's HMP**. The model therefore predicts **0 % xylose-derived MP2P in every pot and no
   MP2P at all in the no-thiamine pot F**, where Table 6 shows it formed at >95 % xylose. The
   "thiamine-diagnostic" label holds only at C&B's dilution (77–90 % thiamine). A sugar + H2S
   route to 3-mercapto-2-pentanone is missing. The paper names the candidate acceptors:
   1-deoxypentosone, norfuraneol, 5-hydroxy-2,3-pentanedione and 2,3-pentanedione. This is the
   pair to the MFT test the repo already runs, and it fails structurally, not numerically.
2. **The S-methyl of 2-methyl-3-(methylthio)furan is 70 % xylose carbon in the ternary and 64 % in
   the no-thiamine pot.** **[!] STRUCTURAL CONTRADICTION.** The model makes `MMFT` only as
   MFT + `MESH`, and `MESH` only from `r_thi_mesh` (thiamine). It therefore predicts an
   all-thiamine methyl in the ternary and **no MMFT in pot F**, where it is formed. A
   sugar-derived C1 (formaldehyde; cf. the 13 % singly labelled thiamine bicycle, Table 4) plus
   cysteine sulfur is the methylthio source when cysteine is present. HOLD-OUT, qualitative.
3. **The B50 thiamine → H2S step at pH 5, 145 °C, the bracket.** In the **no-cysteine pot D, FFT
   forms and is >99 % xylose-labelled**. Sulfur therefore reached a xylose skeleton from thiamine
   at pH 5. B50's share is **0 below pH 6.25** (from Dwivedi & Arnold 1971 at ~100 °C), and pot D
   has no other sulfide source. **The model predicts FFT = 0 in pot D, which is a presence
   contradiction.** The other side of the bracket comes from the same pot: **MFT, 3-mercapto-
   2-pentanone and 2-methylthiophene are each <1 % xylose-labelled**. So whatever sulfide
   thiamine releases is "by far much too small to contribute significantly to the formation of
   MFT" (author). This gives a HOLD-OUT pair for any pH-5 extension of B50 at 145 °C: **FFT > 0
   AND (sugar-route MFT / total MFT) < 0.01 in pot D**. The existing
   `cerny2007b_control_no_cysteine` test checks only the second half, which B50's zero passes
   trivially. It cannot be converted into a share number: there are no amounts, and the sulfur
   may move by S-transfer from a thiamine fragment rather than as free H2S. This points the same
   way as Dwivedi & Arnold 1971 at 120 °C, pH 6.0, where H2S was already 22.5 mmol/mol, not 0.
4. **The concentration → route-mix trend** (§5) with Wang 2012 as a third point. This is a sign
   test of the first-order vs second-order structure. The C&B pair is already a HOLD-OUT.
5. **The thiamine-only markers** (4,5-dihydro-2-methyl-3-furanthiol, the dihydrofuranone, the
   thiophenone, 3-mercapto-2-butanone: ≥95–99 % thiamine) say that the sugar route does not feed
   them. The model carries none of them, so they are nothing to test now.

**Cannot anchor:** any concentration, yield, rate, Ea, pH or time dependence (single pH,
temperature and time, no amounts); H2S; the relative MFT amounts between pots (not reported, so
"cysteine increases MFT" cannot be read here). 2-Methylthiophene's 46 % unlabelled share without
thiamine (a cysteine-carbon route, or an artefact of the M−1 ion basis) has no species in the
model.

---

## 7. CAVEATS

1. **Isotope ratios only.** Nothing in the paper is an amount, so a share of a tiny pool and a
   share of a large one look the same.
2. **Controls are single, unreplicated experiments**; the ternary's replication is not stated.
   The "<1", "<5", ">95" and ">99" entries are detection-limited bounds, not values.
3. **98 % enrichment, and only M and M+5 read** (M+4 for 3-mercapto-2-butanone, M−1-based for
   2-methylthiophene). Scrambled isotopomers (M+1…M+4 from fragment recombination) are not in
   Tables 3/5/6, but Table 2 shows they are small for MFT, MP and FFT (m/z 115–117 / 119–121 ≤
   3.6 %).
4. **Thiamine at 100 mM** is far above food levels (meat ~0.003–0.03 mM [Z], from 0.1–1 mg/100 g), and Cerny himself warns
   the share is concentration-dependent (§5). The 46 % is a property of this loading, which is
   why the repo correctly made the pair, not the single number, the test.
5. **pH 5.00 is the cold buffer pH.** 0.5 M phosphate should hold it, but the pH at 145 °C and
   after the reaction is not measured. Thiamine HCl (100 mM) and cysteine add acid load to the
   buffer.
6. **Xylose, not ribose**; the model's PENT lumps both.
7. Compounds 8 and 12 are tentative identifications.

# Cui, Ma, Wang, Hayat, Zhang & Ho 2022 (JAFC 70:8409) — EXTRACTION
### *"Temperature-Dependent Catalysis of Glycylglycine on Its Amadori Compound Degradation to Deoxyosone"*
### Purified xylose-glycylglycine Amadori compound (XGG-ARP) alone (20 mM) or with added glycylglycine (GG) at ARP : GG = 2 : 1, 1 : 1, 1 : 2, 1 : 3 **at a constant TOTAL of 20 mM** (so the ARP charge falls 13.3 -> 5 mM as GG rises), water, initial pH 7.0 (25 °C), 100 and 120 °C, 0-120 min: XGG-ARP, the cross-linking product XGG-CP (HPLC-ELSD, purified standards), 3-deoxyxylosone (OPD, semi-quantitative, ug/mL) and furfural (SPME). **The "catalysis" claim rests on absolute 3-DX levels in pots whose ARP charge differs 4x; 1-DX, GO, MGO were not observed.**

**Read 2026-10-31 from `data/articles/cui2022.pdf`** (8 pp., born-digital ACS; stamped "Downloaded via SEOUL NATL UNIV
on November 9, 2022"). Prose and Table 1 from `pdftotext -layout`; Figs. 1a/b, 2, 4, 5a/b/c read from 300 dpi renders
cropped per panel (read error ~ +-0.3 mM on ARP, +-0.01 mM on XGG-CP, +-1 ug/mL on Fig. 2, +-0.7 ug/mL on Fig. 5b,
+-0.15 units on Fig. 5c). Supporting Information (FTIR spectra only, Figs. S1-S2) not on disk; it holds no data the
model needs.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Heping Cui and Mengyu Ma (equal contribution), Ziyan Wang, Khizar Hayat, Xiaoming Zhang* (Jiangnan University), Chi-Tang Ho* (Rutgers) |
| Venue | *J. Agric. Food Chem.* **2022**, 70 (27), 8409-8416 |
| DOI | **10.1021/acs.jafc.2c03427** (printed; Crossref 2026-10-31: same title, six authors, vol. 70, issue 27, pp. 8409-8416, issued 2022-06-30) |
| Dates | received 16 May 2022; revised 20 Jun 2022; accepted 21 Jun 2022; published 30 Jun 2022 |
| Laboratory | **Jiangnan University (Zhang X / Cui H) with Ho (Rutgers)** -- the Ho network |
| Companion | Ma et al. 2021 JAFC 69:14915 (ref 18: XGG-ARP / XGG-CP preparation and structure) -- not on disk |

## 1. What this paper answers for the model

| model need | what Cui 2022 gives | tag |
|---|---|---|
| (a) 1-DX / 3-DX split, flux or level | **Neither.** Only 3-DX was observed: "only 3-deoxyxylosone was observed" [V]; no LOD for 1-DX. Standing levels, ex situ OPD (method "similar to" refs 20-21, not described) | [Q] / [!] |
| (b) deoxyosone lifetimes | Not measured. Lower bounds on the 3-DX exit from its decline where the ARP is nearly gone: 120 °C, ARP : GG 1 : 3 pot, 60 -> 120 min **k >= 0.027 /min (lifetime <= ~37 min)**; ARP-alone pot k >= 0.012 /min (<= ~87 min) | [Z] |
| (c) amine (dipeptide) catalysis of ARP -> deoxyosone | **Confounded.** Per ARP CHARGED, peak 3-DX rises with the GG share at BOTH temperatures (100 °C: 1.9 % alone -> 4.3 % at 1 : 3; 120 °C: 0.72 % -> 2.7 %), and the ARP's fractional loss in 30 min rises with GG at 120 °C (k 0.058 -> 0.089 /min) but not at 100 °C (0.026-0.028, one outlier 0.038). The ARP charge (5-20 mM) and GG's own buffering (amine pKa ~8.2, my value) change with the GG share, so neither "catalysis" nor "inhibition" is isolated | [F] / [Z] |
| (d) base-favoured 3-DX sink / amine trap of 3-DX | **None shown.** GG does not lower 3-DX per ARP; the authors' XGG-CP route makes 3-DX, it does not consume it | [F] / [P] |
| (e) mechanism of the "catalysis" | XGG-CP (ARP + GG cross-link) as an intermediate to 3-DX [P]. **XGG-CP never exceeds 0.30 mM (1.5 % of the ARP)**, so it cannot carry a material 3-DX flux unless it turns over very fast (not shown) | [F] / [Z] |
| (f) observation operator | OPD derivatisation, presumably ex situ (protocol not printed); 3,4-hexanedione bought (internal standard implied); **no 3-DX standard**; 3-DX in ug/mL -> semi-quantitative (CF-062 class) | [M] / [!] |

**Answers, weighed.**
1. **The paper does not show amine catalysis of enolisation** (confidence ~0.75). Its two observations -- GG
   "inhibits" 3-DX at 100 °C and "catalyses" it at 120 °C -- are both statements about absolute levels in pots where the
   ARP charge is 20, 13.3, 10, 6.7 or 5 mM. Normalised to the ARP charged, GG raises the 3-DX yield at both temperatures
   by 2-4x [Z]; the temperature dependence the title claims disappears on that normalisation. Rival explanations of the
   per-ARP rise that the design cannot exclude: (i) a 3-DX exit that grows with [ARP] (Zhou 2024b's ARP + 3-DX
   condensation) is weaker at low ARP charge; (ii) less acid made at low ARP charge plus GG buffering keeps the pH higher
   (pH not measured); (iii) genuine dipeptide catalysis.
2. **3-DX is short-lived at 120 °C** (lifetime tens of minutes, an upper bound of ~37-87 min [Z]), and is only ~1 % of
   the ARP charged at its peak: a quasi-steady intermediate, as hub H1 needs.
3. **"Dehydration ... dominant fate of 3-deoxyxylosone" [V] is contradicted by the paper's own Table 1** [Z][!]:
   furfural at 100 °C, 120 min, ARP alone is 5.35 ng/mL = 0.056 umol/L, ~3 x 10^-6 of the 18 mM ARP consumed and
   ~1/4000 of the peak 3-DX.

## 2. Conditions [M]

| item | as printed |
|---|---|
| XGG-ARP | Xyl : GG 2 : 1 (51 g/L), initial pH 8.0, 70 °C 15 min (simultaneous vacuum dehydration, ref 18); Dowex 50WX4 (H+); purity **97.37 %** |
| XGG-CP | Xyl : GG 3 : 1 (55 g/L), pH 7.5, 90 °C 15 min with vacuum dehydration, redissolved, 120 °C 60 min in a pressure bottle; purity **92.15 %** |
| Model pots | ARP alone 20 mmol/L; ARP : GG 2 : 1, 1 : 1, 1 : 2, 1 : 3 "to a total concentration of 20 mmol/L" [V] (Fig. 1a confirms ARP at t = 0 of 13.3 / 10 / 6.7 / 5 mM). Figures label them as GG : ARP 1 : 2, 1 : 1, 2 : 1, 3 : 1. **pH 7.0 at 25 °C, no buffer**, not re-measured. 20 mL split in two: 100 °C and 120 °C; sampled 30-120 min. Vessel and headspace for the model pots not stated |
| ARP, CP | RP-HPLC-ELSD, XBridge BEH amide; log-log curves with the purified products: ARP log y = 0.6889 log x - 4.6523 (R2 0.999), CP log y = 0.4932 log x - 3.6631 (R2 0.992) |
| Dicarbonyls | HPLC-DAD; conditions "similar to those described previously" [V] (refs 20-21, Jiangnan); timing of OPD addition, reagent, IS not printed [!] |
| Furfural | SPME CAR/PDMS, 3 mL + 2 mL sat. NaCl, 50 °C 40 min; external standard y = 1.689e5 x + 7e6 (R2 0.9968; units of x not printed) |
| Replicates | triplicate; ANOVA / Duncan |

## 3. Quantitative results

**ARP (mmol/L) [F]:**

| min | alone 100 | 2:1 100 | 1:1 100 | 1:2 100 | 1:3 100 | alone 120 | 2:1 120 | 1:1 120 | 1:2 120 | 1:3 120 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 20 | 13.3 | 10 | 6.7 | 5 | 20 | 13.3 | 10 | 6.7 | 5 |
| 30 | 9.0 | 5.75 | 4.6 | 2.15 | 2.15 | 3.55 | 1.8 | 0.75 | 0.55 | 0.35 |
| 60 | 5.35 | 3.35 | 2.45 | 2.3 | 1.07 | 2.45 | ~0.35 | ~0.35 | 1.45 [!] | ~0.3 |
| 90 | 3.45 | 2.05 | 2.05 | 1.65 | 0.6 | ~0.75 | ~0.3 | ~0.2 | ~0.75 | ~0.1 |
| 120 | 2.0 | 1.3 | 0.57 | 1.3 | 0.37 | ~0.5 | ~0.25 | ~0.15 | ~0.25 | ~0.1 |

(Columns are ARP : GG. Values below ~0.5 mM at 120 °C overlap within the read error. The ARP : GG 1 : 2 pot at 120 °C
rises 0.55 -> 1.45 mM from 30 to 60 min [!].) [Z] 30-min fraction left: 100 °C 0.45 / 0.43 / 0.46 / 0.32 / 0.43;
120 °C 0.18 / 0.135 / 0.075 / 0.08 / 0.07. Two-point apparent Ea of ARP-alone loss (0-30 min) ~47 kJ/mol [Z].

**XGG-CP (mmol/L) [F]:** 100 °C (Fig. 1b) alone 0.193 / 0.24 / 0.287 / 0.30; 2:1 0.123 / 0.152 / 0.198 / 0.202;
1:1 0.064 / 0.120 / 0.148 / 0.152; 1:2 0.065 / 0.083 / 0.088 / 0.084; 1:3 0.040 / 0.061 / 0.067 / 0.064.
120 °C (Fig. 4) alone 0.25 / 0.27 / 0.28 / 0.26; 2:1 0.17 / 0.185 / 0.178 / 0.157; 1:1 0.13 / 0.136 / 0.125 / 0.11;
1:2 0.077 / 0.066 / 0.060 / 0.058; 1:3 0.048 / 0.041 / 0.036 / 0.034 (30 / 60 / 90 / 120 min).

**3-DX (ug/mL; / 132.12 for mmol/L) [F]:**

| min | alone 100 | 2:1 100 | 1:1 100 | 1:2 100 | 1:3 100 | alone 120 | 2:1 120 | 1:1 120 | 1:2 120 | 1:3 120 |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 50.2 | 47.3 | 26.0 | 27.0 | 28.1 | ~19.1 | 16.1 | ~19.3 | 21.7 | 17.6 |
| 60 | 43.4 | 35.2 | 22.9 | 26.5 | 20.2 | 10.1 | 10.7 | 9.5 | 8.3 | 17.4 |
| 90 | 20.1 | 27.6 | 27.6 | 24.3 | 17.6 | 8.8 | ~9.3 | 6.9 | 14.4 | 9.7 |
| 120 | 15.4 | 27.3 | 20.3 | 20.1 | 13.3 | 5.05 | ~8.1 | 7.9 | ~8.1 | 3.4 |

(At 120 °C, 30 min, the ARP-alone and 1 : 1 points overlap at ~19.2; assignment +-0.7.) [Z] peak 3-DX as % of ARP
charged: 100 °C 1.9 / 2.7 / 2.0 / 3.1 / 4.3 %; 120 °C 0.72 / 0.91 / 1.46 / 2.5 / 2.7 % (alone, 2:1, 1:1, 1:2, 1:3).

**Table 1, furfural at 100 °C (ng/mL) [M]:** alone 2.89 / 3.84 / 4.13 / 5.35; GG : ARP 1:2 0.80 / 1.32 / 2.15 / 3.42;
1:1 0.64 / 1.24 / 1.69 / 1.65; 2:1 0.56 / 1.20 / 1.55 / 2.74; 3:1 0.34 / 0.66 / 1.35 / 1.57 (30-120 min).

**Fig. 5c, furfural at 120 °C, axis printed "ug/uL" [F][!]:** alone 2.7 / 2.9 / 4.3 / 4.9; GG : ARP 1:2 1.4 / 2.25 /
4.2 / 5.55; 1:1 2.2 / 2.9 / 5.15 / 6.9; 2:1 1.45 / 2.55 / 4.9 / 7.95; 3:1 5.0 / 5.7 / 6.85 / 10.5. **At face value
4.9 ug/uL = 51 mM from a 20 mM ARP -- impossible**; in Table 1's ng/mL it would be no higher than at 100 °C. Only
within-figure ratios are usable (3 : 1 / alone = 2.1 at 120 min). The legend markers of 5c differ from the other panels
(1 : 1 square, 3 : 1 triangle); the series above follow 5c's own legend.

**Not observed [Q]:** 1-DX, GO, MGO, pyrazines.

## 4. Mechanism evidence and its kind

| claim | kind | strength |
|---|---|---|
| XGG-ARP + GG -> XGG-CP | purified, NMR / MS-identified product (ref 18), calibrated time course | moderate (identity from the companion paper) |
| XGG-CP -> 3-DX | [P]; no fed-CP experiment in this paper | none |
| GG catalyses ARP -> 3-DX at 120 °C | within-study factor, **confounded** with ARP charge and pH | weak |
| 3-DX -> furfural is the dominant 3-DX fate | [P]; contradicted by the furfural / 3-DX molar ratio | none |
| ARP undergoes 1,2-enolisation only (no 1-DX) | negative observation without LOD | weak |

**Laboratory:** Jiangnan University (Zhang X, Cui H) with Ho (Rutgers).

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| Evidence map 2.1, amine acceleration of enolisation: NONE | **LEAVES UNCHANGED** | the GG factor is confounded with ARP charge and buffering; the hub review's lead ("glycylglycine catalyses ARP -> deoxyosone at 120 C", hub.md sec. 4 item 3) is not supported as stated |
| Hub H3 (split shifts to 1-DX with pH) | LEAVES UNCHANGED | one pH; no 1-DX |
| Hub H1 (engine 3-DX too long-lived) | **SUPPORTS** | 3-DX lifetime <= ~37-87 min at 120 °C [Z], vs the engine's ~35 h at 100 °C |
| B7 base-favoured 3-DX sink (prereg sec. 19-22) | LEAVES UNCHANGED | no pH series, no trap, no amine-dependent 3-DX loss per ARP |
| Peptide amines raise DX 2-8x (evidence map 2.1, Ho-network row) | **QUALIFIES** | the rise exists per ARP charged, but the design cannot attribute it to the amine |
| CF-062 (deoxyosones semi-quantitative) | **SUPPORTS** | no 3-DX standard; units ug/mL |
| CF-060 / CF-071 / CF-076 | LEAVE UNCHANGED | different ARP (dipeptide), no cysteine, no TTCA |

## 6. Disagreements found

1. **Fig. 5c units vs Table 1** (furfural, same pot design, 100 vs 120 °C): ng/mL vs "ug/uL", a factor 10^6 at face
   value; Fig. 5c's face value exceeds the ARP charge 2.5x. Candidate causes: axis-label error (ng/uL or ug/L); a
   different calibration; not settled by the paper.
2. **1-DX absent here, present in every other Jiangnan ARP pot** (Zhou 2022 / 2024b AX-ARP: 1-DX 0.2-0.8 mM; Hou 2026 and
   Wei 2022 GX-ARP: 1-DX 1-15 % of ARP). Candidates: a different ARP (dipeptide N), detection limit (3-DX itself peaks at
   only ~0.15-0.38 mM here), or the HPLC gradient.
3. **3-DX yield per ARP an order of magnitude below the AX-ARP pots of the same laboratory at 120 °C** (peak 0.72 % of
   ARP alone vs 10 % in Zhou 2024b and 30 % in Zhou 2022) -- different ARP, so chemistry is possible, but the CF-060
   method spread is the same order.
4. **Self-contradiction on 3-DX fate** (sec. 1, item 3) [!].

## 7. Ledger-ready measurements (two readers)

1. Figs. 1a / 5a XGG-ARP at 100 and 120 °C, five charges (calibrated with the purified ARP): ARP lifetime vs charge and
   GG share -- the one clean dataset here; tests whether ARP loss is first order in ARP.
2. Figs. 2 / 5b 3-DX (semi_quantitative, shape only): the 3-DX rise-and-fall at two temperatures; per-ARP normalisation
   needs the ARP charge column.
3. Figs. 1b / 4 XGG-CP (calibrated): bounds on the cross-link pool (<= 1.5 % of ARP).
4. Table 1 furfural (calibrated, 100 °C) -- only with the CF-style flag on Fig. 5c's unit.

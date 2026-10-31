# Wei, Cui, Hayat, Zhang & Ho 2022 — EXTRACTION
### *"Glycine-Xylose Amadori Compound Formation Tracing through Maillard Browning Inhibition by 2-Threityl-thiazolidine-4-carboxylic Acid Formation from Deoxyosone and Exogenous Cysteine"*
### Purified glycine-xylose Amadori compound (GX-ARP, 5 mM) +- cysteine (5 mM), water, initial pH 7.0, 110 / 120 / 130 °C, 0-100 min: GX-ARP decay (second-order fits, Ea), 1-DX / 3-DX / GO / MGO (OPD quinoxalines, 120 °C only), cysteine and TTCA (UPLC-MS). Plus a two-step xylose (0.2 M) + glycine (0.1 M) browning experiment (80 °C then 120 °C with 5 mM cysteine). **ONE pH only: signature (vii)'s pH dependence is untested here.** The TTCA-from-deoxyosone pathway is proposed, not shown, and is not mass-balanced.

**Read 2026-10-31 from `data/articles/wei2022.pdf`** (8 pp., JAFC 70:12164-12171; born-digital, clean text layer;
Figs. 1-3 and Scheme 2 checked against 150 / 300 dpi page renders; figure values read by me, read error ~ +-0.5
percentage points on Fig. 3a / 3b, ~ +-0.05 mM on Figs. 2 / 3c).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Question (signature (vii), `results/validation/kinetic_core_dpo_partition_prereg.md` sec. 14; hub hypothesis H2 and
finding 2 in `docs/guides/MECHANISM_REVIEWS_2026-10-31.md` sec. 1): does free cysteine trap the deoxyosones, which one
faster (3-DX vs 1-DX), how does that depend on pH, at what rate?

| need | what Wei 2022 gives | tag |
|---|---|---|
| (a) pH dependence of trapping | **None.** All model pots at initial pH 7.0 (+-0.1), unbuffered water, pH not held or re-measured | [M] |
| (b) which deoxyosone falls more with cysteine (120 °C, 5 mM ARP + 5 mM Cys vs ARP alone; mol % of ARP) | **1-DX falls first:** at 10 min 2.6 vs 5.0 (-48 %), 30 min 6.3 vs 13.9 (-55 %). **3-DX falls later:** 10 min 6.0 vs 6.0 (0 %), 20 min 8.5 vs 13.7 (-38 %), 30 min 15.4 vs 20.7 (-26 %), then declines to 7.7 vs 32.7 at 90 min (-76 %). At 90 min both are -76 %; in absolute terms 3-DX loses more (25 vs 11.5 points) | [F] |
| (c) a trapping rate constant | **None printed.** Rough [Z] only: TTCA rises ~0.0048 mM/min (20-100 min, 120 °C) against ~0.8-1.4 mM free Cys and ~0.4 mM 3-DX + ~0.23 mM 1-DX -> k ~ 0.005-0.01 mM^-1 min^-1 (~5-10 M^-1 min^-1) if ALL TTCA came from DX; an upper bound given sec. (e) | [Z] |
| (d) temperature | TTCA % of ARP at 100 min: **5.0 / 9.6 / 22.4** at 110 / 120 / 130 °C; DX only at 120 °C | [F] |
| (e) is the product of the trap TTCA? | **Not shown.** TTCA identified by UPLC-Q-TOF (Fig. S3, not on disk). Scheme 2 makes TTCA from 1-DX and 3-DX, but a deoxyosone (C5H8O4) + Cys - H2O gives C8H13NO5S (M 235), TTCA is C8H15NO5S (M 237): **the scheme needs a net 2-H reduction it does not supply** (the drawing turns 3-DX's C-3 CH2 into CHOH) | [P] / [Z] / [!] |
| (f) ARP decay with vs without Cys | second-order fits 1/c = 1/c0 + kt, c0 = 5 mM, k in L mmol^-1 min^-1: without Cys **0.01302 / 0.01748 / 0.02403**, with Cys **0.01153 / 0.0148 / 0.02228** at 110 / 120 / 130 °C; Ea 39.32 vs 44.75 kJ/mol | [M] |
| (g) cysteine loss | free Cys 5 -> **~1.55 / 1.42 / 1.24 mM in 20 min** at 110 / 120 / 130 °C, ~0.85 / 0.76 / 0.66 mM at 100 min: ~70 % gone in 20 min at all three temperatures | [F] |

**Answer for signature (vii): no test.** Wei 2022 has one pH (7.0), so it cannot show a trap that "grows with pH".
Confidence the paper holds no pH series: ~99 %. What it adds, weighed:

1. **Cysteine lowers both deoxyosones at pH 7, 120 °C; 1-DX responds first, 3-DX later and more in absolute terms**
   [F]. This mildly contradicts the premise of H2 that the aldehydic 3-DX is the preferred target: at 10 min 3-DX is
   untouched while 1-DX is halved. Weight: **low (~0.3 that it bears on H2 at all)**, because (i) the DX levels are
   net of formation and decay, and cysteine also slows the ARP decay (k -15 % at 120 °C), so a lower level does not
   isolate trapping; (ii) ~70 % of cysteine is gone by 20 min (by routes the paper does not identify), so most of the
   3-DX drop (30-90 min) happens at ~1 mM free Cys; (iii) one experiment series per condition (n = 3 replicates,
   SD bars), one temperature for the DX panels.
2. **The TTCA here most likely does not come from the deoxyosones** (my reading, **~0.65**): Scheme 2 is 2 H short
   (sec. 1 (e)). A mass-balanced route exists: Amadori reversal to the xylose Schiff base / xylose (shown to be a
   major Amadori fate at 37 °C, zyzak1995_extraction.md: hexoses ~20 % of fructoselysine at half-life) followed by
   transimination / condensation with cysteine, or Cys-ARP <-> TTCA interconversion (zhai2020_extraction.md). The
   authors exclude free xylose only because "xylose was not observed" [Q] (no detection limit given; free xylose would
   be scavenged by Cys or Gly). Consequence for the engine: Wei's TTCA rows, if keyed, measure ARP -> (reversal) ->
   TTCA, not DX trapping. They would then support an ARP-reversal source of TTCA in the Zhai pots, ~10 % (120 °C) to
   ~22 % (130 °C) of 5 mM ARP in 100 min [F].
3. **The cysteine sink is mostly not TTCA** [Z]: at 120 °C, 100 min, 4.24 mM Cys consumed vs 0.48 mM TTCA (11 %); at
   130 °C 4.34 vs 1.12 mM (26 %). The fast, nearly temperature-independent 70 % drop in 20 min looks like oxidation
   (air, no deaeration stated, no chelator in the reaction pots) or adduct formation, not an activated Maillard step.
   This bears on the cysteine-fate rival "headspace O2 oxidises thiol" (prereg sec. 12, P ~0.2): a dilute-pot
   cysteine loss of ~70 % in 20 min, flat across 110-130 °C, is what that rival predicts. Weak support (~0.55).
4. **ARP decay numbers** (usable as a check on the engine's GX-ARP lifetime, water, pH 7 initial) [M]: at c0 = 5 mM the
   initial pseudo-first-order rate k c0 is 0.065 / 0.087 / 0.120 min^-1 (t1/2 = 15 / 11 / 8 min) at 110 / 120 /
   130 °C without Cys [Z]. The second-order form is an empirical fit (R^2 0.946-0.997), not a mechanism. Recomputing
   Ea from the printed k: 39.3 kJ/mol without Cys (matches 39.32), **42.2 kJ/mol with Cys (paper prints 44.75)** [Z][!].
5. **Deoxyosone yields from GX-ARP alone (pH 7, 120 °C)** [F]: 1-DX peaks ~14 % of ARP at 30-40 min, 3-DX rises to
   ~33 % at 90 min; GO ~12.6 % and MGO ~5 % (peak) -- a 3-DX / 1-DX ratio of ~1.5 at 30 min and ~2.2 at 90 min in water
   at pH 7. These are direct data for the hub's DX partition (H1 / H3), at one pH.

## 2. Identity

| field | value as printed |
|---|---|
| Authors | Shangjie Wei, Heping Cui* (Jiangnan University), Khizar Hayat, Xiaoming Zhang, Chi-Tang Ho* (Rutgers) |
| Title | as in the header |
| Venue | *J. Agric. Food Chem.* **2022**, 70 (38), 12164-12171 |
| DOI | **10.1021/acs.jafc.2c04961** (printed on every page; Crossref 2026-10-31: same title, authors, vol. 70, issue 38, pp. 12164-12171) |
| Dates | received 13 Jul 2022; revised 7 Sep 2022; accepted 8 Sep 2022; published 20 Sep 2022 |
| Content | Figs. 1 (A420 vs cysteine addition time; GX-ARP at 80 °C), 2 (A420; GX-ARP decay -Cys / +Cys at 110-130 °C), 3 (1-DX, 3-DX, Cys, TTCA, GO, MGO), Schemes 1-2. No tables |

## 3. Methods [M]

| experiment | as printed |
|---|---|
| Two-step browning | xylose 0.2 M + glycine 0.1 M, pH 7.0, deionised water, 80 °C water bath for 20-160 min; ice; + cysteine 5 mM, pH re-adjusted to 7.0 +- 0.1; oil bath 120 °C, 120 min; A420 |
| GX-ARP | xylose 0.2 M + glycine 0.1 M, pH 7.0, 80 °C 90 min; vacuum 60 °C; Dowex 50WX8 (H+), water then 0.1 M NH4OH; purity "99.05%" (Chemicals) / "over 99.5%" (Methods) [!]; NMR (Fig. S1) shows furanose forms (structures Fig. S2) |
| Model pots | GX-ARP 5 mM + cysteine 5 mM (control: GX-ARP alone), pH adjusted to 7.0, **no buffer named** (water) [!]; oil bath 120 °C, sampled 20, 40, 60, 80, 100 min; kinetics at 110, 120, 130 °C. Fig. 3a/b/e/f plot 10, 20, 30, 40, 50, 60, 90 min [!]. Vessel, headspace, deaeration: not stated |
| Dicarbonyls | OPD derivatisation (reagent: HEPES 1 M, "DTPA (11 mol/L)" [!: impossible concentration, a typo], OPD 2 %), 12 h dark RT; HPLC-DAD 315 nm, C18; method of Yu et al. (ref 21) |
| TTCA, Cys, xylose, ARP | UPLC-ESI-Q-TOF (positive), BEH C18; quantification method / standards not described [!] |
| Kinetics | 1/c = 1/c0 + kt fitted in GraphPad Prism; Arrhenius over 110-130 °C; triplicates, ANOVA / Duncan |

## 4. Results

**Fig. 1 [F]:** A420 after the two steps is lowest (~0.48) when cysteine is added at 80 min (0.67-0.76 at 50-70 min,
0.61-0.67 at 90-100 min); GX-ARP at 80 °C peaks at ~15.3 mM at 80 min (from 0.1 M glycine: ~15 %), then ~12.5 mM.

**Fig. 2b / 2c, GX-ARP (mM) [F]:**

| min | 110 -Cys | 120 -Cys | 130 -Cys | 110 +Cys | 120 +Cys | 130 +Cys |
|---|---|---|---|---|---|---|
| 0 | 5 | 5 | 5 | 5 | 5 | 5 |
| 20 | ~2.05 | ~1.78 | ~1.58 | ~1.97 | ~1.73 | ~1.53 |
| 40 | ~1.38 | ~1.10 | ~0.88 | ~1.50 | ~1.32 | ~0.95 |
| 60 | ~1.07 | ~0.82 | ~0.60 | ~1.32 | ~1.07 | ~0.68 |
| 80 | ~0.92 | ~0.65 | ~0.42 | ~1.13 | ~0.88 | ~0.50 |
| 100 | ~0.72 | ~0.52 | ~0.25 | ~1.08 | ~0.78 | ~0.42 |

Printed fits (k, L mmol^-1 min^-1; R^2): -Cys 0.01302 (0.9968), 0.01748 (0.9464), 0.02403 (0.9778); +Cys 0.01153
(0.9620), 0.0148 (0.9871), 0.02228 (0.9986) [M].

**Fig. 3, 120 °C, mol % of initial GX-ARP [F]** (panel identities from the y-axis labels and the text; see flag on
the caption):

| min | 1-DX -Cys | 1-DX +Cys | 3-DX -Cys | 3-DX +Cys | GO -Cys | GO +Cys | MGO -Cys | MGO +Cys |
|---|---|---|---|---|---|---|---|---|
| 10 | 5.0 | 2.6 | 6.0 | 6.0 | 3.0 | 3.1 | 1.1 | 0.7 |
| 20 | 10.2 | 4.6 | 13.7 | 8.5 | 6.9 | 6.3 | 1.9 | 0.8 |
| 30 | 13.9 | 6.3 | 20.7 | 15.4 | 8.1 | 8.0 | 2.9 | 1.4 |
| 40 | 14.0 | 5.6 | 22.9 | 12.2 | 9.1 | 3.0 | 4.1 | 0.7 |
| 50 | 13.4 | 5.0 | 24.6 | 8.7 | 10.3 | 0.9 | 5.1 | 0.7 |
| 60 | 12.6 | 4.7 | 25.4 | 8.1 | 10.4 | 3.7 | 4.8 | 1.0 |
| 90 | 15.2 | 3.7 | 32.7 | 7.7 | 12.6 | 6.1 | 4.3 | 0.9 |

**Fig. 3c, free cysteine (mM) [F]:** 110 °C 5 / 1.55 / 1.22 / 1.04 / 0.92 / 0.85; 120 °C 5 / 1.42 / 1.12 / 0.87 /
0.78 / 0.76; 130 °C 5 / 1.24 / 1.05 / 0.74 / 0.69 / 0.66 at 0 / 20 / 40 / 60 / 80 / 100 min.

**Fig. 3d, TTCA (% of GX-ARP) [F]:** 110 °C 1.6 / 2.2 / 2.8 / 4.2 / 5.0; 120 °C 2.0 / 4.8 / 6.3 / 7.6 / 9.6; 130 °C
3.1 / 9.5 / 16.4 / 21.5 / 22.4 at 20 / 40 / 60 / 80 / 100 min.

**Shown vs proposed.** Shown: cysteine lowers A420, DX, GO, MGO and slightly slows ARP loss; TTCA appears and grows
with T; no GX-ARP-derived thiazolidine (Hou 2017's adduct) detected [Q]. Proposed [P]: cysteine's amine condenses with
a DX carbonyl, rearranges, and the thiol closes the ring to TTCA (Scheme 2); GO / MGO also fall by Strecker
degradation with cysteine [P].

## 5. Supplementary material

**Exists, not on disk.** The SI paragraph lists: Fig. S1 (1H / 13C NMR of the purified GX-ARP), Fig. S2 (proposed
isomer structures of the GX-ARP), Fig. S3 (UPLC-Q-TOF spectrum of TTCA). **No time courses, rate constants or tables.**
Only Fig. S3 matters to the model: it would show whether the product is m/z 238 ([M+H]+ of TTCA) or 236 (the
deoxyosone adduct Scheme 2 actually implies), i.e. it can settle flag (e). Worth fetching only for that check.

## 6. Flags

- [!] **No pH series** -- signature (vii) and H2's pH dependence are not addressed.
- [!] **Scheme 2 is not redox-balanced** (deoxyosone + Cys cannot give TTCA without a reduction); TTCA from DX is
  asserted, not traced (no labelling, no fed-DX control).
- [!] **Fig. 3 caption mislabels panels:** caption "3-DX (a), 1-DX (b), MGO (e), GO (f)" and "cysteine (d) and TTCA
  (d)"; axis labels and text say (a) 1-DX, (b) 3-DX, (c) cysteine, (d) TTCA, (e) GO, (f) MGO. This dossier follows
  the axes and text.
- [!] Sampling times in Methods (20-100 min) differ from Fig. 3a/b/e/f (10-90 min).
- [!] Ea with Cys: 44.75 printed vs 42.2 recomputed from the printed k (two-point and three-point fits agree on 42.2).
- [!] Unbuffered water, pH set at start only; headspace / O2 unstated; GX-ARP purity printed two ways (99.05 / 99.5 %).
- [!] TTCA quantification standard and response factor not described.

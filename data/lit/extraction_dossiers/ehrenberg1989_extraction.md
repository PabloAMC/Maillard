# Ehrenberg, Harms-Ringdahl, Fedorcsák & Granath 1989 — EXTRACTION (Cu- and Fe-catalysed oxidation of cysteine by O2: rate law, constants, Fe–Cu interaction)

**Source on disk:** `data/articles/Ehrenberg1989.pdf`, an 11-page image scan (tiff2pdf, 2007; no text layer,
so `pdftotext` returns nothing and no text cross-check was possible). Every number below was read by eye from
the page images, rendered at 250 dpi and cropped table by table (the Table 6 header again at 600 dpi), on
2026-10-09. Journal pages 177–187 are PDF pages 1–11. Written for the question in
`docs/guides/EXPERIMENTS.md` sec. 1: is the engine missing a transition-metal-catalysed oxidation of cysteine?

| field | value |
|---|---|
| Title | "Kinetics of the Copper- and Iron-Catalysed Oxidation of Cysteine by Dioxygen" |
| Authors | L. Ehrenberg, M. Harms-Ringdahl, I. Fedorcsák (Dept. of Radiobiology), F. Granath (Dept. of Mathematical Statistics), University of Stockholm |
| Venue | Acta Chemica Scandinavica 1989, 43, 177–187 (received June 27, 1988; experiments done 1974–1975, p. 186) |
| DOI | 10.3891/acta.chem.scand.43-0177 (not printed on the scan; taken from the citation as supplied) |

## 1. Methods (p. 186, "Experimental")

- L-cysteine HCl, Tris HCl, EDTA disodium, Ellman's reagent (Sigma); CuCl2·2H2O, FeSO4·7H2O,
  FeNH4(SO4)2·12H2O, NaH2PO4, Na2HPO4, H2O2, TiOSO4 (Merck).
- Water: deionised, distilled in glass, redistilled twice in quartz. Glassware washed with **0.05 M EDTA**,
  rinsed in Cu-free distilled water, dried at 200 °C.
- **Cleanliness criterion**: 1 mM cysteine in the buffer under efficient aeration had to stay unchanged
  "within the accuracy of the analysis (± 2 %), over 40 min at 37 °C".
- Most runs at **37.0 °C** (some 29.0 °C) in **40 mM Tris, pH 7.2 or 8.1, measured at the experimental
  temperature**; O2 held at equilibrium with air or with pure O2 (CO2-free, water-saturated) by bubbling with
  vigorous stirring; usually 30 ml. Fe stocks made in 0.05 M H2SO4 and neutralised with NaOH at the start.
  Cu in stocks by atomic absorption; Fe(II) spectrophotometrically at 304 nm after H2O2 oxidation.
- RSH by Ellman's reagent; H2O2 by Ti(IV) (Marklund's modification of Bonet-Maury); O2 uptake by Warburg in
  parallel runs. k2′ measured anaerobically (N2 or Ar).
- k0 and K fitted by least squares on numerical solutions of eqns (6a, b) (eqn 13; MLAB). Time unit: minute.

## 2. Findings that matter

### 2a. The scheme and the rate law (pp. 177–179)

Printed equations (p. 178 eqns 1–4; p. 179 eqns 5–6e; p. 180 eqn 7):

| eqn | as printed |
|---|---|
| (1) | 2 RSH + O2 → RSSR + H2O2 (Cu-catalysed; "only the second of which is spontaneous" refers to (2)) |
| (2) | 2 RSH + H2O2 → RSSR + 2 H2O |
| (3) | 4 RSH + O2 → 2 RSSR + 2 H2O |
| (4) | 2(1+a)RSH + O2 → (1+a)RSSR + (1−a)H2O2 + 2aH2O, 0 < a < 1 |
| (5) | d[RSH]/dt = − k2′ [RSH][H2O2] |
| (6a) | d[RSH]/dt = − k0′ − k2′ [RSH][H2O2] |
| (6b) | d[H2O2]/dt = k0′/2 − (k2′/2)[RSH][H2O2] |
| (6c) | k0′ = k0[RSH]/(K + [RSH]) |
| (6d) | [H2O2]SS = k0′/(k2′[RSH]) |
| (6e) | lim[RSH]→∞ d[RSH]/dt = − 2k0′ |
| (7) | k0[RSH]/(K+[RSH]) = c[RSH]/(1 + [RSH]/K), c = k0/K, "the semi-first order rate constant for reaction (1) at low concentrations of RSH ([RSH] ≪ K; under the conditions studied K ≤ 1×10⁻³ M)" |
| (8) | (Zwart et al., pH ≈ 13.5) k0 (= −d[O2]/dt) = kI[O2]^0.5×[Cu] + kII[O2]^0.5×[Cu]² |
| (9)–(10) | RSH + Cat (+O2) ⇌ RSH–Cat(–O2) (ka, kb) → Cat + Products (kc); K = (kb + kc)/ka |
| (11) | [H2O2]SS = k0/(k2′(K+[RSH])) (→ k0/(k2′[RSH])) |

Abstract (p. 177): the Cu(II)-catalysed reaction (1) "follows Michaelis–Menten kinetics with respect to RSH,
and probably O2, and is, partly at least, second order with respect to Cu". Reaction (2) is first order in
RSH and H2O2, enhanced by Fe(II)/Fe(III) "but practically not at all by Cu(II)". With Fe(II) or Fe(III)
as catalyst the overall reaction proceeds without H2O2 formation and is first order in RS⁻ and in Fe(II)/Fe(III)
"(at least at high [Fe])"; an FeO(OH) sol makes it zeroth order in RSH. "The Cu- and Fe-catalysed rates are
not additive."

**Thiol:Cu complex.** No stoichiometry measured. The Discussion (p. 184) proposes a Cu(II) complex,
"presumably Cu(SR)2", citing EPR showing "100 (±2) % of the Cu to be present as Cu(II)"; stability constants
"it has not been possible to determine". A binuclear (Cu+Cu or Cu+Fe) complex is suggested (pp. 183–185).

### 2b. Reaction (2), cysteine + H2O2 (Table 1, p. 178; text p. 178)

| T / °C | buffer | pH | k2′ / M⁻¹ min⁻¹ | k2 / M⁻¹ min⁻¹ |
|---|---|---|---|---|
| 25 | Tris 40 mM | 7.2 | 160 | 1780 |
| 25 | | 7.7 | 410 | 1810 |
| 25 | | 8.2 | 830 | 1660 |
| 25 | | 8.7 | 1300 | 1710 |
| 25 | | 9.2 | 1640 | 1800 |
| 25 | | | mean | 1730 |
| 37 | Tris 40 mM | 7.2 | 310 | 2780 |
| 37 | | 7.7 | 757 | 2700 |
| 37 | | 8.3 | 1584 | 2590 |
| 37 | | 8.7 | 2090 | 2620 |
| 37 | | | mean | 2700 |
| 37 | Phosphate 10 / 60 / 120 mM | 7.2 | 350 / 380 / 400 | — |

k2′ refers to total cysteine, k2 to the thiolate RS⁻ (Table 1 footnote a: pKa 8.2 at 25 °C and 8.1 at 37 °C).
Text: "k2 ≈ 2700 M⁻¹ min⁻¹ at 37 °C ... k2 ≈ 1700 M⁻¹ min⁻¹ at 25 °C, corresponding to **ΔH* ≈ 27 kJ mol⁻¹**".
Check (derived here): Ea = R ln(2700/1700)/(1/298.15 − 1/310.15) = 29.6 kJ/mol; ΔH‡ = Ea − RT ≈ 29.6 − 2.5 =
27.1 kJ/mol, matching the printed value. **This is the only temperature coefficient in the paper and it
belongs to reaction (2), not to the catalysed reaction (1).**

### 2c. Catalysts and EDTA on reaction (2) (Table 2, p. 179; 40 mM Tris, pH 7.2, 37 °C)

| catalyst | catalyst / µM | k2′ / M⁻¹ min⁻¹ |
|---|---|---|
| none | — | 310 |
| Cu(II) | 0.5–1 | 340 |
| Fe(II) | 0.3 | 735 |
| Fe(II) + 0.3 µM EDTA | (0.3) | 800 |
| Fe(II) + 3 µM EDTA | (0.3) | 310 |
| Fe(II) | 1 | 1 600 |
| Fe(II) | 3 | 5 500 |
| Fe(II) | 10 | 14 000 |

**The EDTA result**: EDTA at 1:1 with Fe(II) did not inhibit (800 vs 735); a ten-fold excess abolished the
catalysis ("This catalytic action of Fe is completely inhibited by a ten-fold excess of EDTA", p. 179). The
paper has no chelator arm on the Cu-catalysed reaction (1).

### 2d. Cu-catalysed autoxidation, 37 °C, 40 mM Tris pH 7.2 (Table 3, p. 180)

| [RSH]0 / mM | [Cu] / µM | interval / mM | −d[RSH]/dt / µM min⁻¹ | d[H2O2]/dt / µM min⁻¹ | k0′ / µM min⁻¹ | [H2O2] mid, expected / µM | found / µM |
|---|---|---|---|---|---|---|---|
| 30 | 1 | 30–20 ᵇ | 240 | 0 | 120 | 14 | (15) ᶜ |
| 10 | 1 | 10–5 | 2.3(1)×10² ᵈ | 2.3 | 117(10) ᵈ | 49 | 48 |
| | | 5–1.7 | 150 | 2.2 | 90 | 59 | 70 |
| 3 | 1 | 3–1.4 | 148 | 5 | 80 | 88 | 92 |
| 1 | 1 | 1–0.5 | 70(5) ᵈ | 15 | 51 | 80 | 70 |
| 1 | 0.5 | 1–0.4 | 32 | 2.5 | 19 | 52 | 50 |
| 1 | **0.25** | 1–0.7 | **14** | 0.6 | **8** | 24 | 18 |
| 0.3 | 1 | 0.3–0.2 | 35(1) ᵈ | 11 | 29 | 46 | 44 ᵉ,ᶠ |

Footnotes: ᵃ intervals within which the rates are linear in time; ᵇ "Solution becomes turbid at about 20 mM
because of precipitation of cystine"; ᶜ from an experiment with cysteamine; ᵈ mean and maximum error of three
experiments; ᵉ rapid rise of [H2O2]; ᶠ final [H2O2] = 104 µM, "70 % of the amount expected to be generated
(150 µM)". Oxygen is not stated in the caption; Fig. 1 (same conditions) says "bubbled with air".

Fig. 2 (p. 179) is a Lineweaver–Burke plot of 1/k0′ vs 1/[RSH] at [Cu(II)] = 1 µM and 0.3 µM, pH 8.1, 37 °C,
in O2 — the basis for the Michaelis form (6c).

### 2e. pH, O2 and [Cu] dependence of k0′ (Table 4, p. 180; 37 °C, 40 mM Tris; k0′ from eqn 6e)

k0′ / µM min⁻¹ at [Cu] / µM =

| PO2 | pH | [RSH] / mM | 0.3 | 1 | 2.5 | 5 | 10 |
|---|---|---|---|---|---|---|---|
| Air | 7.2 | 1 | 12 | 67 | 85 | 100 | 118 |
| | | 5 | 30 | 81 | 163 | 191 | 225 |
| | 8.1 | 1 | 3 | 48 | 210 | 720 | – ᵃ |
| | | 5 | 8 | 60 | 334 | 800 | 870 |
| O2 | 7.2 | 1 | 14 | 107 | 450 | – ᵃ | – ᵃ |
| | | 5 | 40 | 240 | 740 | 850 | 1 400 |
| | 8.1 | 1 | 5 | 37 | ≈200 | ≈900 | – ᵃ |
| | | 5 | 9 | 97 | ≈550 | ≈1 950 | – ᵃ |

ᵃ "Too fast for measurement."

Text (p. 181): at pH 8.1 "k0 (and k0′) being proportional to the square of the [Cu], possibly with saturation
at 10 µM Cu. At pH 7.2 a saturation effect predominates at [Cu] ≥ 1 µM." O2: "Michaelis' constants
**KO2 ≈ 0.05 atm at pH 8.1 and ≈ 1.4 atm at pH 7.2**". Discussion (p. 184): at [Cu] ≈ 1 µM the rate shows a
maximum at pH 7.2–7.4 (refs 7, 8); at pH ≈ 8 "the reaction order with respect to Cu was found to be decidedly
> 1, probably around 2, at [Cu] = (0.3–5) µM". Recalculating Zwart's pH-13.5 data to pH 7.2, 1 µM Cu, air,
37 °C "generates a rate that is some 20 times lower than that observed".

### 2f. k0/K (Table 5, p. 181) — a unit discrepancy

(k0/K) / min⁻¹ (s.e.) at [Cu] / µM = 1, 2.5, 5:

| PO2 | pH | [RSH] / mM | 1 | 2.5 | 5 |
|---|---|---|---|---|---|
| Air | 7.2 | 1 | 21(2) | 43(9) | 39(17) |
| | | 5 | 18(4) | – ᵃ | 24(5) |
| | 8.1 | 1 | 5(0.3) | 32(1) | 116(5) |
| | | 5 | 6(0.3) | 38(1) | 116(8) |
| O2 | 7.2 | 1 | 41(5) | 136(40) | – ᵇ |
| | | 5 | 41(7) | 74(6) | 103(17) |
| | 8.1 | 1 | 7(1) | 36(2) | 157(8) |
| | | 5 | 7(2) | 19(3) ᶜ | 138(18) |

ᵃ not followed to completion; ᵇ too fast; ᶜ "obviously too low by a factor ca. 2". Data at 0.3 µM Cu not
analysable.

**Caution.** Fig. 3 (p. 181) plots "(k0/K)/min⁻¹" for the same means, and the plotted values (read from graph,
approx.: pH 7.2 air ≈ 0.19 / 0.38 / 0.31; pH 8.1 ≈ 0.06 / 0.3 / 1.15 at 1 / 2.5 / 5 µM) are Table 5's means
divided by about 100. Fig. 3's scale is the physically consistent one: from Table 3, k0 ≈ 120 µM min⁻¹ (30 mM
row) and k0′ = 29 µM min⁻¹ at 0.3 mM give K ≈ 0.94 mM and k0/K ≈ 0.13 min⁻¹ (derived here: K = 120×0.3/29 −
0.3), close to Fig. 3, not to Table 5. **Do not use Table 5's numbers in min⁻¹ as printed.**

### 2g. Fe-catalysed autoxidation (Table 6, p. 182; Figs 4–7)

First-order constants for cysteine of initial concentration 1–5 mM, 40 mM Tris, averages of 3–6
determinations over six months, maximum variation about 10 %. **The column header is printed "k2′/min⁻¹"**; the
caption, text and footnote d call these first-order constants k1′, which is what they are.

| added catalyst | [Fe] / µM | PO2 | pH 8.1, 29 °C | pH 8.1, 37 °C | pH 7.2, 37 °C |
|---|---|---|---|---|---|
| none | 0 | O2 | <0.0002 | <0.0002 | <0.0002 |
| Fe(II) | 1 | Air | – | – | ≤0.001 ᵃ |
| | | O2 | 0.009 | – | – |
| | 10 | Air | 0.090 | 0.075 | 0.016 ᵇ |
| | | O2 | 0.17 ᶜ | 0.14 | 0.024 ᵇ |
| Fe(III) | 1 | Air | 0.003 | – | – |
| | | O2 | 0.006 | – | – |
| | 10 | Air | 0.065 ᵈ | – | 0.015 |
| | | O2 | 0.15 | 0.13 | 0.027 (0.007) ᵉ |

ᵃ reaction order could not be determined; ᵇ tendency to lower rates at [RSH] = 5 mM; ᶜ values as low as
0.13 min⁻¹ occasionally; ᵈ rate constant for O2 consumption: k1′ = 0.067 in 40 mM Tris, 0.054 min⁻¹ without
buffer; ᵉ mean and maximum error.

- **Metal-free background**: k1′ < 0.0002 min⁻¹ at pH 7.2 and 8.1, 29 and 37 °C, in O2.
- First order in RSH (log[RSH] vs t linear at ~10 µM Fe, Fig. 4). Fe(III) about as effective as Fe(II). Rate
  rises monotonically with pH 7.2–9.5, "approximately proportional to [RS⁻]" (Fig. 5). Air → O2 roughly
  doubles the rate. Linear in [Fe] above ca. 5 µM, ∝ [Fe]^1.5 below (Fig. 7), possibly from Fe losses to glass.
- Fig. 5 (Fe(III), 29 °C; [Fe] not printed in its caption; the pH 8.1 point matches Table 6's 10 µM rows), read
  from graph, approx.: k1′ ≈ 0.03 / 0.085 / 0.17 / 0.245 / 0.30 min⁻¹ at pH ≈ 7.2 / 7.6 / 8.1 / 8.5 / 9.0; the
  thiolate constant k1 = k1′(1 + 10^(pKa−pH)) is flat at ≈ 0.33–0.35 min⁻¹.
- Fig. 7 (Fe(II)), read from graph, approx.: k1′ rises from ≈ 0.025 min⁻¹ at ≈ 3 µM to ≈ 0.17 at 10 µM and
  ≈ 0.3–0.45 at 20–25 µM.
- No H2O2 detectable; catalase had no effect; "within a few per cent the disappearance of 4 mol equivalents RSH
  was concomitant with the consumption of 1 mol of O2" (p. 185), i.e. eqn (3). Tris did not affect the rate.
- Aged/hydrolysed Fe(III) (FeO(OH) sol) shifts kinetics toward zeroth order in RSH (Fig. 6).
- **Temperature**: the only 29 vs 37 °C pairs (pH 8.1) give 37 °C/29 °C ratios of 0.075/0.090 = 0.83,
  0.14/0.17 = 0.82, 0.13/0.15 = 0.87 (derived here): no positive temperature coefficient over 8 K.

### 2h. Mixed Cu + Fe (Table 7, p. 183; Fig. 8; text pp. 183–185)

40 mM Tris pH 7.2 (one row pH 8.1), equilibrium with air, 37 °C.

| [RSH]0 / mM | [Cu] / µM | k0′ / µM min⁻¹ | −d[RSH]/dt at [Fe] = 0 / 0.1 / 0.3 / 1 / 10 µM (µM min⁻¹) | [H2O2]t½ at [Fe] = 0 / 0.1 / 0.3 ᵈ (µM) |
|---|---|---|---|---|
| 0.3 | 1 | 35 | 42 / – / 37 / 42 / – | 56 / – / 0 |
| 1 ᵇ | 0.2 | 6.3 | 12 / – / 6.5 / 6.2 / – | 40 / – / 0 |
| 1 ᵃ | 0.86 | 45 | 63 / 40 / 34 / 39 / 88 | 130 / 44 / 14 |
| 1 ᵇ | 5 | 102 | 137 / – / 160 / 170 / 500 | 156 / – / 0 |
| 5 | 1 | 90 | 145 / 100 / 63 / 50 / 180 | 85 / – / – |
| 1 (pH 8.1) | 1 | 18 | 34 / 16 / 12 / – / – | 17 / 8 / 0 |
| 1 ᵇ,ᶜ | 1 | 56 | 72 / 34 / 28 / 28 / – | 116 / 32 / 0 |

ᵃ mean of 5 experiments, C.V. < 5 %; ᵇ mean of two; ᶜ FeSO4 stock made in water; ᵈ "At 1–10 µM Fe, no H2O2
could be determined." (The [H2O2] header prints "[Fe]/µm".)

Text (p. 185): "With 1 µM Cu(II) + (0.3–1) µM Fe(II), −d[RSH]/dt assumes a value which is about 40 % of that
obtained with 1 µM Cu(II) alone", kinetics move to stricter zeroth order (smaller K), and H2O2 disappears. At
0.3–1 µM, Fe alone has hardly any catalytic effect on reaction (1) (Fig. 8 caption: "Fe alone at these
concentrations has no influence on [RSH]"). At 10 µM Fe the rate rises again (Table 7). Cu and Fe are "not
independent, i.e. not additive".

### 2i. Matrix effects (pp. 185–186)

In complete media 0.1 and 0.5 mM in each of nineteen other amino acids: at 0.1 mM, k0′ and k2′ "practically
unchanged"; at 0.5 mM "a certain reduction", attributed to Cu(II)–cysteine–histidine complexes competing with
the catalytic Cu–Cys complex. E. coli (10⁷ cells/ml): slight changes.

## 3. What it means for the model

**What the engine runs today for cysteine alone in buffer** (live values; `engine._B2_FIT_REPORT` resolves to
`results/validation/kinetic_core_b9_fit_report.json`, confirmed in the `maillard` env on 2026-10-09):

| step | live value | source |
|---|---|---|
| `r_cys_thermal` (`k_cys_thermal`), first order in Cys, no pH, no metal, no O2 term | log10 k(145 °C) = **−2.066359667088082** /min (k = 8.58×10⁻³ /min); key `b8.k_cys_thermal.log10_k_ref_145C`, centre −2.066359667088082 in `results/validation/core_prediction_uncertainty.json` (the frozen optimum, sampled normal_log10 around it) | b9 `frozen_parameters.log10_k_ref_at_145C` |
| its barrier | **Ea = 55.1 kJ/mol**, fixed by `MEASURED_EA_OVERRIDES` | `KANG_EA_FREE_CYS_DEPLETION_KJ_MOL`, `src/kinetic_core/parameters_sulfur.py:1898`; anchor `KANG_CYS_ANCHOR` (Zhai 2023 / Kang 2026: 10 mM Cys, pH 7, sealed, 100–140 °C) |
| `r_cys_h2s` (`k_cys_h2s`), measured | pH-5 pair A = 1.93×10¹² s⁻¹, Ea = 133.0 kJ/mol | `ZHENG_CYSTEINE_THERMOLYSIS`, `parameters_sulfur.py` (Zheng & Ho 1994) |
| `ch_cys_ox` (`k_cys_ox`), Cys + dissolved O2 | **0 (inert)**: b9 carries no `oxygen` block, so `engine.shipped_oxygen_consumers()` returns 0.0. The uncertainty file's `sulfur.oxygen.k_cys_ox.log10_k` centre 0.0 is a placeholder for "not a free coordinate" (`uncertainty.py` ~l. 988), **not** log10 k = 0. Even if active it has no metal term | `src/kinetic_core/sulfur.py:777` |
| `k_thiolate_loss` | acts on FFT and MFT only, not on Cys | `sulfur.py:687–702` |

Derived here: k_cys_thermal(95 °C) = 8.58×10⁻³ × exp[−(55.1/R)(1/368.15 − 1/418.15)] = 9.97×10⁻⁴ /min;
k_cys_h2s(pH 5, 95 °C) = 1.93×10¹² × 60 × exp(−133.0/(R·368.15)) = 1.56×10⁻⁵ /min; sum 1.013×10⁻³ /min;
survival exp(−5 × 1.013×10⁻³) = **0.9949 at 5 min** and exp(−180 × 1.013×10⁻³) = **0.833 at 3 h**. This
reproduces EXPERIMENTS.md's 99.5 % and 83 %, so these are the constants behind that probe. At 37 °C the same
two steps give 3.44×10⁻⁵ /min (1.1 µM/min at 33 mM).

**What Ehrenberg predicts at the EXPERIMENTS.md conditions, as far as the paper reaches.** Measured
conditions: 37 °C, pH 7.2, 40 mM Tris, air bubbled continuously. Target: 33 mM Cys, 0.25 µM Cu.

1. At high [RSH] the Cu channel is **zero order in cysteine** (eqn 6e: −d[RSH]/dt → 2k0′, K ≲ 1 mM) and
   **catalyst-limited**. Table 3 gives k0′ = 8 µM/min at 1 mM / 0.25 µM Cu and 51 at 1 mM / 1 µM; scaling the
   saturated 30 mM / 1 µM value (120) by 8/51 gives k0′ ≈ 18.8 µM/min, total ≈ 2k0′ ≈ 38 µM/min (derived here).
   The measured 1 mM / 0.25 µM total, **14 µM/min**, is the lower edge. Band: **14–38 µM/min**.
2. At 33 mM that is 70–190 µM in 5 min, **0.2–0.6 %** of the charge (derived here: 5 × 14 / 33 000 and
   5 × 37.6 / 33 000); zero-order half-life 7.3–20 h; effective pseudo-first-order 4.2×10⁻⁴–1.1×10⁻³ /min. That is
   12–33× the engine's 37 °C rate (3.44×10⁻⁵), and about equal to what the engine already does at 95 °C
   (1.0×10⁻³ /min).
3. Losing "nearly all" (say 90 %) of 33 mM in 5 min needs ≈ 5.9 mM/min, 160–420× the 37 °C, pH 7.2 rate. From
   37 °C that factor needs Ea ≈ 190–230 kJ/mol by 60 °C, or 83–99 kJ/mol by 95 °C (derived here:
   Ea = R ln(f)/(1/310.15 − 1/T)). Because the reaction is zero order, the fraction lost scales as 1/[Cys]:
   Table 3's own 0.3 mM / 1 µM Cu row (35 µM/min) clears that pot in about 10 min at 37 °C. **"Within minutes"
   is what this rate law gives at sub-millimolar cysteine and ~1 µM Cu, not at 33 mM and 0.25 µM.** The
   cysteine concentration of the "real pot" EXPERIMENTS.md cites is not given in that paragraph and should be
   checked against its source.
4. **Oxygen caps the channel in a sealed vial.** Eqn (1) uses 1 O2 per 2 RSH; eqn (3) 1 per 4 (Fe measured
   ≈ 4:1, p. 185). In EXPERIMENTS.md's 20 mL vial with 5 mL liquid: cysteine 165 µmol; dissolved O2 at about
   0.25 mM (a general value, not from this paper) is ≈ 1.25 µmol, enough for 1–3 % of the cysteine; the 15 mL of
   air headspace holds ≈ 128 µmol O2 (derived here: 0.2095 × 15×10⁻⁶ m³ × 101 325 Pa / (R × 298.15 K)), enough
   for 257–514 µmol RSH. So the channel needs gas–liquid transfer, and every rate in this paper was measured
   with bubbling. At pH 7.2 KO2 ≈ 1.4 atm, so the rate is close to proportional to pO2 at air.
5. **Cystine solubility**: Table 3 footnote b reports turbidity from cystine at about 20 mM. A 33 mM pot that
   oxidises appreciably will precipitate cystine.

**Extrapolation to 95–145 °C and to pH 5 needs a declared Ea band and a declared pH law.** The paper offers
only the following:

- **Temperature**: no Ea or ΔH‡ for the catalysed reaction (1). ΔH* ≈ 27 kJ/mol is for reaction (2) alone. The
  Fe channel shows no increase from 29 to 37 °C (ratios 0.82–0.87). Any Ea for the Cu channel would be a
  declared assumption, and it has to be combined with O2 solubility, which falls with temperature.
- **pH**: data at pH 7.2 and 8.1 only. At ~1 µM Cu the rate peaks at pH 7.2–7.4 (refs 7, 8) and falls at 8.1. The
  Fe channel goes with the thiolate fraction (Fig. 5). If the Cu channel also went with the thiolate fraction
  (pKa 8.1 at 37 °C), pH 5 would scale the 7.2 rate by
  [1/(1+10^(8.1−5))]/[1/(1+10^(8.1−7.2))] = 7.9×10⁻⁴/0.112 = **7.1×10⁻³** (derived here). At saturation
  (zero order, Michaelis-bound complex) it could fall less. Kreitman 2016 (see `kreitman2016_extraction.md`)
  shows that Cu(II) is still reduced to a Cu(I)–thiolate instantly at pH 3.6, so complexation survives low pH;
  turnover does not. A pH factor between ~10⁻² and 1 relative to pH 7.2 is the honest range from these papers.
- **Fe at pH 5** (derived here, an extrapolation below the measured pH 7.2): k1′ ≈ 0.34 min⁻¹ ×
  1/(1+10^3.15) ≈ 2.4×10⁻⁴ min⁻¹ at ~10 µM Fe, 29 °C.
- **Chelator**: EDTA 1:1 with Fe does nothing, 10:1 abolishes (Table 2; reaction 2 only). Amino acids at 0.5 mM
  measurably slow the Cu channel (p. 186), and a Maillard pot is full of them.
- **Fe + Cu**: 0.1–1 µM Fe cuts the 1 µM Cu rate to about 40 % and removes H2O2; 10 µM Fe raises it again.

## What it does not give

- No Ea or ΔH‡ for the Cu- or Fe-catalysed autoxidation; nothing above 37 °C; nothing below pH 7.2.
- No measured thiol:Cu complex stoichiometry or stability constant (the paper says it could not determine them).
- No phosphate-buffer data for the catalysed reaction (phosphate only in Table 1, for reaction 2), and no chelator
  arm on the Cu reaction. Table 5's k0/K units are inconsistent with Fig. 3 by about 100× (sec. 2f).
- No measurement in a sealed vessel; all O2 supply was by bubbling.

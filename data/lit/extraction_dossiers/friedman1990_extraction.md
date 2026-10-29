# Friedman & Molnar-Perl 1990 — EXTRACTION (β-alanine 400 + glucose 200 mmol/L, 100 °C, 15-120 min, pH ~6, closed 5 mL vials; browning as A420 "molar absorptivity"; L-cysteine, N-acetyl-L-cysteine, GSH, bisulfite, urea as inhibitors; NO H2S, NO sulfide, NO cysteine or sugar loss, NO rates)

### NOT A SOURCE FOR THE STRECKER H2S SHARE. The paper measures only brown colour (A420, plus spectra at 200-600 nm and a little NMR). Sulfur volatiles, H2S, sulfide, thiol consumption, cysteine loss and glucose loss are never measured or mentioned. The US 5,760,005 citation rests on the closing sentence about parenteral amino-acid + carbohydrate solutions, which is an aspiration, not data.

**Source on disk:** `data/articles/friedman1990.pdf` (6 pp., journal pp. 1642-1647; scanned, Acrobat Capture OCR).
Read with `pdftotext -layout -q`. Pages 1643-1644 (Tables I-III, Fig. 2) re-read from 200-dpi page renders through the
Read tool, because the Table II text layer is scrambled (columns run together; "0.750.25" etc.). The numbers below come from
the images.

## 0. Identity

| field | value |
|---|---|
| Title | "Inhibition of Browning by Sulfur Amino Acids. 1. Heated Amino Acid-Glucose Systems" ✔ |
| Authors | Mendel Friedman, Ibolya Molnar-Perl — Western Regional Research Center, USDA-ARS, Albany CA (Molnar-Perl visiting from L. Eötvös University, Budapest) ✔ |
| Venue | J. Agric. Food Chem. **38**(8), 1642-1647 (1990). Received 12 Jun 1989, revised 21 Mar 1990, accepted 3 Apr 1990 ✔ |
| DOI as printed | **None printed** in the PDF (no DOI string anywhere in the text layer or on the page images). Not filled in here. |
| Companion papers | Parts 2 (fruit juices, protein foods) and 3 (apples, potatoes), Molnar-Perl & Friedman 1990a,b, same issue. Not on disk as far as this read goes |
| Paper type | Browning-inhibition screen: end-point A420 at a fixed time, with an "index of prevention". Single-wavelength colour only |
| Cited by | US 5,760,005 (as the source on heated amino-acid + glucose solutions with cysteine / NAC) |

## 1. Why it matters (and why it mostly does not)

It was downloaded in case it measured H2S or sulfide from cysteine (free amine) vs N-acetylcysteine (blocked amine)
with glucose. That would have constrained wave B59a's release of cysteine's Strecker sulfur as H2S, whose floor comes from Kobayasi 1965
(4.2-8 mol % of cysteine at 100 °C, 1 h). **It does not.** What it delivers:

1. A **browning (A420) time course and concentration series** for β-alanine + glucose at 100 °C (Table II), plus α-alanine
   and Nα-acetyl-L-lysine end points (Table I).
2. A **pH profile of browning** at 100 °C, 90 min, initial pH 4.4-10.0 (Fig. 2, graph only).
3. **Minimum thiol doses for ~70-90 % browning suppression** at 100 °C, 90 min (Table III): L-cysteine works at 0.05 mol/mol,
   NAC needs 0.2 mol/mol. So the free-amine thiol is about 4× more potent per mole than the N-blocked one. That is the
   only cysteine-vs-NAC contrast in the paper, and it is a colour observable.
4. A **null result**: N2 vs O2 saturation gave the same browning, with or without NAC (no numbers given).

## 2. Methods as they matter to a kinetic model

- **Stocks** (all adjusted to pH 6.5 before use): D-glucose 2 M; β-alanine 2 M; Nα-acetyl-L-lysine 0.5 M; amino-acid mix 0.6 M;
  glycylglycine 0.8 M; NAC 0.5 M; L-cysteine (free base) 0.4 M; GSH 0.4 M; NaHSO3 0.1 M; urea 5 M.
- **Pots**: 0.5 mL glucose stock + 1.0-2.0 mL amine stock + inhibitor, made to **5.0 mL with distilled water**, in **screw-cap
  vials** (Kimble 60910-1), so the vessel was closed with a small headspace. **No buffer** (pH set on the stocks only).
- **Standard inhibitor conditions**: 100 °C (boiling-water bath), **90 min**, pH 5.8-6.1. Glucose **200 mmol/L**;
  β-alanine **400 mmol/L**; Nα-acetyllysine **100 mmol/L**; glycylglycine **160 mmol/L**; amino-acid mix **120 mmol/L**.
  Amino-acid mix (mol/L in its stock): α-Ala 0.02; Arg, Gly, His, Ile, Orn, Phe, Pro, Ser, Thr, Val 0.04 each.
  The paper also writes "0.06 mol/L mixed amino acids" for it. The listed components sum to 0.42 mol/L, so the two figures do not agree (flag 3).
- **pH drift**: Table I gives initial pH 6.1-6.2 and final 5.7-5.8 (β-Ala, 90 min). Fig. 2: initial 4.40-10.0 maps to final 4.40-6.85.
  Initial ≈ final up to pH ~7; above ~8 the final pH falls a lot.
- **Analytics**: A420 against water blanks in 1 cm cells (Beckman DB). Reported as "molar absorptivity" Am/(M·cm), but
  **the normalising concentration is not stated** (flag 1). Standard error from preliminary runs is about ±5 %.
- **Never measured**: H2S, sulfide, any volatile, residual cysteine or NAC, thiol titre, glucose, amine loss, intermediates
  (3-DG etc.). No temperature other than 100 °C. No rate constants, no Ea.

## 3. Numbers, re-typed (verbatim from the page images)

### Table I — Am/(M·cm), 100 °C, 90 min, initial pH 6.1-6.2, final 5.7-5.8; amine/Glu = 0.4 M/0.2 M (β-AL), 0.10 M/0.2 M (N-ALL)

| reactants | 420 nm | 300 nm | Am300/Am420 |
|---|---|---|---|
| DL-Al–D-Glu | 0.28 ± 0.015 | 7.0 | 25 |
| L-Al–D-Glu | 0.28 ± 0.015 | 7.0 | 25 |
| D-Al/D-Glu | 0.46 ± 0.010 | 11.0 | 24 |
| β-Al–D-Glu | 3.50 ± 0.20 | 52.5 | 15 |
| N-ALL–D-Glu | 1.43 ± 0.20 | 41.0 | 29 |

(The footnote gives only the β-AL and N-ALL concentrations. The α-alanines are presumably 0.4/0.2 M, as Fig. 1 states "AL/D-Glu = 0.4 M/0.2 M".)

### Table II — Am/(M·cm) at 420 nm, 100 °C, conditions "see Table I"

β-AL:D-Glu time course

| mole ratio | concn, mol:mol (as printed; read as mol/L) | 15 min | 30 min | 60 min | 120 min |
|---|---|---|---|---|---|
| 1:1 | 0.5:0.5 | 0.28 | 1.63 | 0.75 | |
| 1:2 | 0.5:1.0 | 0.47 | 1.90 | 2.25 | |
| 1:3 | 0.5:1.5 | 0.81 | 3.70 | 5.00 | |
| 2:1 | 1.5:0.5 | 1.02 | 3.60 | 5.00 | |
| 3:1 | 1.5:0.5 | 3.00 | 9.50 | 11.50 | |
| 3:1 | 0.75:0.25 | 0.88 | 3.40 | 5.00 | 22.80 |
| | 0.375:0.125 | 0.16 | 0.50 | 0.96 | 5.60 |
| | 0.150:0.050 | 0.10 | 0.16 | 0.24 | 0.90 |
| | 0.075:0.025 | 0.0 | 0.0 | 0.0 | 0.60 |
| 3:1 | 0.60:0.020 | | | 4.40 | |
| 6:1 | 0.60:0.10 | | | 3.40 | |
| 12:1 | 0.60:0.05 | | | 3.00 | |
| 24:1 | 0.60:0.025 | | | 2.92 | |

Amine:glucose series at 90 min (glucose 0.2 M)

| mole ratio | concn, mol:mol | β-AL:D-Glu (90 min) | N-ALL:D-Glu (90 min) |
|---|---|---|---|
| 0.25:1 | 0.06:0.2 | 0.55 ± 0.05 | 1.08 ± 0.08 |
| 0.5:1 | 0.1:0.2 | 0.61 ± 0.04 | 1.43 ± 0.02 |
| 1:1 | 0.2:0.2 | 0.90 ± 0.03 | 3.25 |
| 1.5:1 | 0.3:0.2 | 2.07 ± 0.05 | 5.28 ± 0.20 |
| 2:1 | 0.4:0.2 | 3.50 ± 0.20 | 8.50 |
| 3:1 | 0.6:0.2 | 7.93 ± 0.35 | |
| 4:1 | 0.8:0.2 | 15.30 ± 0.40 | |
| 6:1 | 1.2:0.2 | 18.50 ± 0.10 | |

### Table III — inhibitors; reactions I β-Ala, II Nα-acetyllysine, III amino-acid mix, IV glycylglycine (with glucose; 100 °C, 90 min)

| reaction | inhibitor | IP, % | stoichiometry of inhibition, mol inhibitor / mol MP |
|---|---|---|---|
| I | N-acetyl-L-cysteine | 70 | 0.2 |
| II | N-acetyl-L-cysteine | 83 | 0.4 |
| III | N-acetyl-L-cysteine | 91 | 2.0 |
| IV | N-acetyl-L-cysteine | 89 | 0.4 |
| I | L-cysteine | 79 | 0.05 |
| II | glutathione | 83 | 0.08 |
| I | sodium bisulfite | 79 | 0.02 |
| II | sodium bisulfite | 96 | 0.12 |
| III | sodium bisulfite | 74 | 0.16 |
| IV | sodium bisulfite | 91 | 0.05 |
| I-IV | urea | 91 / 88 / 95 / 89 | 12 / 25 / 12 / 12 |

IP = 100 − MAV(amine + glucose + inhibitor) × 100 / MAV(amine + glucose). "MP, Maillard reaction product precursor."
In the text (abstract and p. 1644) the minimum ratios are given **per mole of D-glucose**: bisulfite 0.02, L-cysteine 0.05,
NAC 0.2, GSH 0.2, urea 8. These do not match the table's GSH (0.08) or urea (12-25) entries (flag 2).
L-cysteine was tested in reaction I only.

### Fig. 2 — pH profile (graph only; values read by eye, ± ~10 %, broken nonlinear y-axis)

β-AL/Glu 0.4/0.2 M, 100 °C, 90 min: Am at 420 nm ≈ 1 at initial pH 4.4-4.9, ≈ 3-4 at pH 6.1, ≈ 10 at about pH 8.5, and ≈ 40-80 at pH 9.4-10.0.
The curve is sigmoid. These are graph reads, **not** tabulated values.

### Text-only statements
- Inhibitor ratios of 0.2, 0.5 and 0.6 (sulfur amino acids or sulfite to MP) gave A420 that **rose with inhibitor dose**.
  The authors attribute this to inhibitor autoxidation products absorbing.
- NMR (β-Ala ± inhibitors): no spectral change with cysteine, NAC, GSH or bisulfite. Urea forms an adduct with β-Ala.
- N2 vs O2 saturation (β-Ala + Glu ± NAC): "the extent of browning of the two solutions was the same".

## 4. What the model could take — graded

Grades: **M** = measured, usable as printed; **F** = fit/derived by the authors; **Z** = my arithmetic (shown); **!** = flag / do-not-use.

| # | item | grade | value / note |
|---|---|---|---|
| 1 | H2S / sulfide from cysteine vs NAC with glucose | **!** | **Not measured.** Not a measurement of the Strecker H2S share, not a bound on it. Unrelated to B59a |
| 2 | cysteine / NAC / glucose loss | **!** | Not measured |
| 3 | rate constants, Ea | **!** | None reported. 100 °C only, so no Ea can be derived |
| 4 | β-Ala + Glu A420 time course, 0.75:0.25 M, 100 °C | **M** | 0.88 / 3.40 / 5.00 / 22.80 at 15 / 30 / 60 / 120 min (Am units, normaliser unknown) |
| 5 | shape of that time course | **Z** | 30→60 min ×1.47, 60→120 min ×4.56. It accelerates (lag/autocatalytic), not first order. A model that browns with a 1st-order intermediate would undershoot at 120 min |
| 6 | concentration dependence (3:1, 60 min) | **Z** | 0.75:0.25 vs 0.375:0.125 → 5.00/0.96 = 5.2 for a 2× dilution, apparent order log2(5.2) = **2.4** in Am. If Am = A/[Glu]·1 cm, absolute A ratio = 10.4, order **3.4**. At 120 min: 22.80/5.60 = 4.1 → 2.0 (or 3.0). Strongly super-linear; the normaliser ambiguity (flag 1) moves it by 1 |
| 7 | amine dependence at fixed Glu 0.2 M, 90 min | **M/Z** | β-Ala 0.06→1.2 M: 0.55→18.50 (×34 for ×20 amine). Between 0.2 and 0.8 M: 15.30/0.90 = 17 for ×4, order ≈ **2.0** in amine (Z) |
| 8 | β-Ala vs α-Ala vs Nα-acetyllysine, 90 min | **M** | 3.50 / 0.28 / 1.43 (β-Ala at 0.4 M; N-ALL at only 0.1 M) |
| 9 | pH sensitivity of browning, 100 °C | **M (graph)** | Am ≈ ×3-4 from pH 4.4 to 6.1, ≈ ×10 more by pH 9.4-10. Qualitative pH check only; initial pH, unbuffered, with drift above pH 8 |
| 10 | cysteine browning suppression | **M** | IP 79 % at 0.05 mol Cys / mol (β-Ala system, Glu 200 mmol/L → Cys **10 mmol/L**), 100 °C, 90 min |
| 11 | NAC browning suppression | **M** | IP 70 % at 0.2 mol/mol (→ NAC **40 mmol/L**), same system |
| 12 | residual A420 with inhibitor | **Z** | Cys: (1 − 0.79) × 3.50 = **0.74**; NAC: (1 − 0.70) × 3.50 = **1.05** Am (β-Ala reference 3.50 from Table I) |
| 13 | free-amine vs blocked thiol potency | **Z** | Cys reaches a higher IP at 1/4 the dose of NAC → ≥ 4× more potent per mol. This is consistent with glucose–cysteine thiazolidine trapping (needs the free NH2) plus SH, but the authors do not test or claim this; NMR saw "no interaction" |
| 14 | O2 independence of browning | **M (qual.)** | N2 = O2; no numbers |

## 5. Flags

1. **"Molar absorptivity" is undefined.** A420 is divided by some molarity that is never named (glucose? amine? total?). Every Am value
   carries this ambiguity. Absolute A420 cannot be recovered with certainty, and apparent reaction orders shift by 1 depending on the choice.
2. **Table III vs text disagree** on the stoichiometric basis (per mol MP vs per mol glucose) and on GSH (0.08 vs 0.2) and urea (12-25 vs 8).
3. **Amino-acid mix concentration inconsistent**: the listed components sum to 0.42 mol/L, yet the text says "0.06 mol/L".
4. **Table II internal oddities** (as printed, image-verified): the 1:1, 0.5:0.5 row **falls** from 1.63 (30 min) to 0.75 (60 min). The "2:1" and "3:1"
   rows both print 1.5:0.5. And "3:1 0.60:0.020" is a 30:1 ratio (probably 0.060:0.020 or a mislabel). Treat those rows as typos; use the
   0.75:0.25 dilution series and the 90-min amine series.
5. **No buffer, no headspace data**; closed screw-cap vials. pH was set only on the stock solutions.
6. **No DOI printed**; none supplied here.
7. **Patent citation overstates it**: US 5,760,005 leans on the paper for parenteral-type cysteine/NAC + amino-acid + glucose
   solutions. The paper has only browning data, in model systems, and mentions parenteral solutions only in the conclusion.
8. Registry: not checked against `data/keys/papers.yml`. β-alanine and acetyllysine are not cysteine systems.

## Use in the model

- **B59a Strecker H2S share: nothing.** It is not a measurement and not a bound. The Kobayasi 1965 floor (4.2-8 %) stands
  unchallenged and unsupported by this paper.
- **At most a qualitative browning check** for a non-sulfur amine. A420 accelerates with time, depends super-linearly on concentration
  (apparent order ~2-3), and rises steeply with pH at 100 °C. That suits a shape test of a melanoidin/colour lane if the engine has one,
  and only after a decision on flag 1. It is not usable as a fit target. Single temperature, one lab, normaliser unknown.
- **Cysteine as a browning inhibitor** (79 % at 10 mmol/L vs 200 mmol/L glucose; NAC needs 40 mmol/L for 70 %) is a qualitative
  check that a cysteine-trapping (thiazolidine / thiol-carbonyl) sink suppresses colour, and that the free amine matters. Evidence class: level_only, colour only.

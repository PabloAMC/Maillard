# Ghiron, Quack, Mawhinney & Feather 1988 — EXTRACTION
### *"Studies on the Role of 3-Deoxy-D-erythro-glucosulose (3-Deoxyglucosone) in Nonenzymatic Browning. Evidence for Involvement in a Strecker Degradation"*
### Synthesised 3-DG (Madson & Feather) with L-alanine / L-phenylalanine at 100 °C: browning rate (A450, pH 3.5 phthalate, 0.1 M), disappearance of 3-DG and alanine (0.4 M each, unbuffered water, 0-7 h, GLC), 13C-NMR of 3-DG + [1-13C] / [2-13C]alanine (no adduct seen), and phenylacetaldehyde IDENTIFIED (not quantified) from 3-DG + Phe. **Second laboratory (Missouri, Feather) for the EXISTENCE of 3-DG as a Strecker donor; no yield. Its Fig. 2 is the only on-disk fed-3-DG + amino acid pot with the DONOR followed over time.**

**Read 2026-10-31 from `data/articles/ghiron1988.pdf`** (4 pp., JAFC 36:677-680; scanned with an OCR layer that
garbles "D-" as "~-" and the figure axes; p. 4 is mostly the following paper by Hsieh et al.). Prose from the text
layer, checked against 200 dpi renders (scratchpad `gh-*.png`); Figs. 1 and 2 read by me from the renders: read
error ~ +-0.03 A450 (Fig. 1), ~ +-3 percentage points and ~ +-0.1 h (Fig. 2).

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 1. What this paper answers for the model

Questions: the forced-input Strecker layer (prereg secs. 1-3) declares **kappa_TDG** (3-DG consumption by the amino
acid) from Martins 2005's 3-DG + Gly melanoidin constant and lists gap **G7**: "a fed GO / 3-DG / glucosone + amino
acid pot with the DONOR followed over time". Does 3-DG give a Strecker aldehyde at all, how much, and how fast is it
consumed by an amino acid?

| model need | what Ghiron 1988 gives | tag |
|---|---|---|
| 3-DG is a Strecker donor | **yes, qualitatively:** PA from 3-DG + Phe, 100 °C, 30 min, identified by GLC co-injection and MS of the oxime-TMS derivative | [Q] |
| Strecker yield from 3-DG | **none printed** | [M] (absence) |
| 3-DG consumption with / without amino acid | **Fig. 2**, 0.4 M each, water, 100 °C, 0-7 h | [F] |
| amino-acid consumption by 3-DG | Fig. 2: alanine loss with 3-DG vs with glucose | [F] |
| a covalent 3-DG-amine intermediate pool | **none detected** by 13C-NMR (D2O, 100 °C, hourly to 5 h) | [Q] |
| pH | browning at pH 3.5 (phthalate); the loss and Strecker pots unbuffered, pH not stated | [M] |
| oxygen / vessel | sealed 6-mm glass tubes; headspace not stated | [M] |

**Answers.**
1. **3-DG is a Strecker donor in a second laboratory** (Missouri vs DFA Garching vs Quest) [Q]. Existence only.
2. **A fed-3-DG pot with the donor followed (gap G7, partially).** From Fig. 2 [F] -> [Z]: 3-DG alone loses 33 %
   in 1 h, 3-DG + Ala 72 %; the alanine-dependent excess first-order rate over 0-1 h is 0.0145 min^-1 at a mean
   [Ala] ~0.36 mol/L, i.e. **kappa_TDG(Ala) ~ 4e-5 L mmol^-1 min^-1** (over 1-2 h: ~1e-5). This is **1.3-1.9 dex
   BELOW** the prereg's declared kappa_TDG (Martins `k_tdg_mel` 8.12e-4 at 100 °C, pH 6.8, Gly) and below the
   prereg's hard floor for TDG (>= 7e-5, from the Glc/Phe pot at pH 7). The pH of Ghiron's unbuffered pots is not
   given (0.4 M Ala + 3-DG in water: plausibly 5-6 and falling [Z]); under the prereg's free-amine law a pH of ~5.5
   instead of 6.8 would account for ~1.2 dex [Z]. So it **neither confirms nor contradicts** the declared value; it
   is a pH-confounded lower anchor. Weight low (~0.3).
3. **Amino acid consumed per extra 3-DG consumed** (the prereg's nu, declared 1, band [0.1, 1]) [Z]: ~0.54 over the
   first hour, ~2.4 by 7 h (alanine keeps falling after the 3-DG is ~93 % gone: 3-DG products also consume it).
   A time-dependent ratio; the declared constant nu is a simplification.
4. **No accumulating Schiff base / adduct** (13C-NMR of [1-13C] and [2-13C]Ala with 3-DG, 5 h at 100 °C) [Q]:
   consistent with no long-lived 3-DG-amine pool at the NMR's detection limit (not printed; a few % [Z]). Weakly
   bears on the CF-073 explanation "a reversible aldehyde-amine pool" (that one concerns PA-amine, not 3-DG-amine;
   unchanged).

## 2. Conditions [M]

| experiment | as printed |
|---|---|
| 3-DG | synthesised (Madson & Feather 1981); purity not stated; analysed as oxime-TMS (mannitol IS), two GLC peaks |
| Browning (Fig. 1) | 0.1 M L-Ala and 0.1 M 3-DG, HMF or glucose, each "made up in 0.1 M phthalate buffer (pH 3.5)"; 0.5 mL + 0.5 mL in a flame-sealed 6-mm tube; boiling water bath; A450 after 10-fold dilution. [!] Mixing equal volumes halves each to 0.05 M, but the caption says "initial concentrations ... all 0.1 M" |
| Loss of reactants (Fig. 2) | 1-mL samples, reactants **all 0.4 M**, **distilled water** (no buffer), sealed 6-mm tubes, 100 °C; aliquot dried under N2 at < 25 °C, derivatised (oxime-TMS for sugars / 3-DG; t-BDMS for Ala, glycerol IS), GLC |
| NMR | 120 mg 3-DG (0.74 M [Z]) + 40 mg Ala + 10 mg [1-13C]- or [2-13C]Ala (~0.56 M total [Z]) in 1.0 mL D2O, dioxane IS, 5-mm tube, 100 °C water bath, spectra hourly to 5 h (~2500 scans) |
| PA identification | 15 mg 3-DG + 14 mg L-Phe in 2.0 mL water (46 mM 3-DG, 42 mM Phe [Z]), sealed, 100 °C, 0.5 h; 3 x 2 mL CHCl3; oxime-TMS; GLC (two peaks = authentic PA oxime-TMS) and GLC-MS (m/z 207, 192, 116, 103, 91, 89, 77, 65) |
| Replicates | not stated; no error bars |

## 3. Quantitative results (all figure-only)

**Fig. 1 — A450 (10-fold diluted), L-Ala + carbonyl, pH 3.5, 100 °C** [F]:

| h | 3-DG | HMF | glucose |
|---:|---:|---:|---:|
| 1 | 0.08 | 0.06 | ~0.01 |
| 2 | 0.31 | 0.08 | ~0.01 |
| 3 | 0.53 | 0.10 | ~0.01 |
| 4 | 0.77 | 0.12 | ~0.02 |
| 5 | 0.96 | 0.13 | ~0.02 |
| 6 | 1.03 | 0.14 | ~0.02 |
| 7 | 1.04 | 0.14 | ~0.02 |
| 8 | 1.04 | 0.17 | ~0.03 |

Sigmoid for 3-DG (lag ~1 h, plateau from ~6 h). "Neither D-glucose nor L-alanine alone gives rise to measurable
color" [Q].

**Fig. 2 — per cent loss of reactant, 0.4 M each, water, 100 °C** [F]:

| h | 3-DG alone | 3-DG + Ala | Ala + 3-DG | glucose + Ala | Ala + glucose |
|---:|---:|---:|---:|---:|---:|
| ~0.15-0.25 | 26-30 | 27 | 7 | ~0 | 0 |
| ~0.4 | — | 49 | 10 | — | — |
| 1 | 33 | 72 | 21 | ~1 | 0 |
| 2 | 46 | 81 | 26 | ~1 | 0 |
| ~4.5 | 68 | 87 | 36 | ~2 | 0 |
| 7 | 74-77 | 93 | 46 | ~4 | 0 |

Derived [Z] (sec. 1): first-order k(3-DG alone) 0.0067 min^-1 (0-1 h), 0.0036 (1-2 h); with Ala 0.0212 / 0.0065;
kappa_TDG(Ala) ~4e-5 (0-1 h) / ~1e-5 (1-2 h) L mmol^-1 min^-1; Ala per excess 3-DG 0.54 (1 h), 2.4 (7 h).
[!] The ~30 % drop of 3-DG ALONE within the first ~15 min, followed by a slower decline, is not first-order; it may
be a fast equilibration (3-DG forms; cf. weenen1992: five major isomers in water) or an artefact of the dry-down /
derivatisation. The early points carry more weight in the 0-1 h kappa than they deserve; the 1-2 h value is the
safer one.

## 4. Mechanism evidence and its kind

Laboratory: **University of Missouri, Columbia (Feather group; Mawhinney GC methods).**
- **Fed intermediate (kind 2):** synthesised 3-DG fed with Ala / Phe; PA identified [Q]; 3-DG and Ala losses
  measured (Fig. 2).
- **Within-study factor (kind 6):** carbonyl identity (3-DG vs glucose vs HMF) on browning and on Ala loss.
- **Isotope-labelled NMR (kind 1 as a negative):** 13C-Ala shows no new C-1 / C-2 signals: no accumulating Schiff
  base or stable adduct [Q]. The Strecker conclusion is by elimination plus PA identification [P for "major
  degradative reaction"]: no CO2, NH3 or aminoketone balance, no PA yield.
- Authors' claim, not shown: the 3-DG + Phe reaction's "major degradative reaction involves a Strecker
  degradation" (abstract). Without a PA amount this is [P].

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| Forced-input Strecker layer SINGLE-LAB (EX-049 / EX-051) | **LEAVES UNCHANGED** for magnitudes; **SUPPORTS** (second lab) only the existence of the 3-DG -> PA step | no yield |
| Declared kappa_TDG (Martins, 8.1e-4) and its hard floor (7e-5) | **LEAVES UNCHANGED (pH-confounded)**; a lower anchor ~1e-5-4e-5 at unknown pH, Ala | if Ghiron's pot was near pH 7 it would contradict both by >= 1.3 dex; its pH is unknown, so carry it as a flagged comparator, not a band edge |
| nu (amine per donor consumed, declared 1, band [0.1, 1]) | **weak CONTEXT** | 0.5 early, > 1 late: the ratio is not constant when donor products also consume amine |
| CF-073 (MGO / 3-DG donor ratio 35 vs 1.6) | **LEAVES UNCHANGED** | no MGO arm, no yield |
| Evidence map 2.2 (3-DX / 3-DG sink) | **adds a 3-DG + amine consumption series** (hexose, water, 100 °C) | the amine accelerates 3-DG loss ~3x in the first hour |
| Hub: 3-DG lifetime alone in water at 100 °C | **CONTEXT** | t1/2 ~2 h after the initial drop [Z]; weenen1992 sees no change at 70 °C / 2 h by NMR, consistent within an Arrhenius factor of ~15-30 [Z] |

## 6. Disagreements found

1. **kappa_TDG: Ghiron-derived ~1e-5-4e-5 vs prereg declared 8.12e-4 L mmol^-1 min^-1** (1.3-1.9 dex) and vs the
   prereg hard floor 7e-5 (0.2-0.8 dex). Candidate causes: pH (unbuffered, unknown; the prereg's free-amine law
   gives ~1.2 dex for pH 5.5 vs 6.8), amine (Ala vs Gly / Phe), concentration (0.4 M: 3-DG hydrate / isomer
   speciation), and the Martins constant being a fitted melanoidin step, not a measured consumption. Not entered as
   a conflict (pH unknown); worth a line in the prereg's kappa_TDG band discussion.
2. **weenen2001_extraction.md cites Ghiron as "3-DG + Phe in water at 100 C, 30 min: Strecker aldehyde formed"**:
   agrees with the PDF (46 mM 3-DG, 42 mM Phe, 0.5 h).
3. Internal [!]: Fig. 1's "all 0.1 M" vs the method's 1:1 mixing of 0.1 M solutions (0.05 M each).

## 7. Ledger-ready measurements

- **Fig. 2** (3-DG alone, 3-DG + Ala, Ala + 3-DG, glucose +- Ala; 0.4 M, water, 100 °C, 0-7 h): the only on-disk fed
  3-DG + amino acid pot with the donor followed; a cross-lab comparator for kappa_TDG and for the hub's 3-DG lifetime.
  Key two readers; carry pH as "unbuffered, not measured".
- **Fig. 1** (A450 sigmoid, 3-DG / HMF / glucose + Ala, pH 3.5, 100 °C): a browning shape, only if the engine ever
  scores colour.
- Nothing for the Strecker yield: PA was identified, not measured.

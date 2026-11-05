# Szilágyi & Györgydeák 1979 — EXTRACTION
DOI: 10.1021/ja00496a026 (not printed in the PDF; Crossref 2026-11-05: "Comments on the putative stereoselectivity in cysteine-aldehyde reactions. Selective C(2) inversion and C(4) epimerization in thiazolidine-4-carboxylic acids", Szilagyi, Gyorgydeak, J. Am. Chem. Soc. 101(2), 427-432, 1979. **Matches the brief and the PDF.**)
### A STEREOCHEMISTRY PAPER, NOT A KINETICS PAPER. (R)- and (S)-cysteine + benzaldehyde or p-tolualdehyde give ~1 : 1 C(2)-epimer mixtures of 2-ARYL thiazolidine-4-carboxylic acids (NMR epimer shares 44-56 % in Me2SO, pyridine, acetic acid, CDCl3, C6D6). N-acetylation inverts C(2) selectively (through a Schiff-base ring opening, proposed); hot acetic anhydride epimerises C(4) WITHOUT ring opening (through a mixed anhydride, proposed). **No rate constant, no equilibrium constant, no water, no sugar ring, no pH.**

**Read 2026-11-05 from `data/articles/szilagyi1979.pdf`** (6 pp., JACS pp. 427-432; scan with a text layer; the last
page continues into an unrelated Fréchet et al. paper). Prose from `pdftotext -layout`; **Table II (solvents and
temperatures, lost in the text layer) read from a 200 dpi render of p. 429** (scratchpad `img/szi-3.png`).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Organic Chemistry and Antibiotics Research Group, Lajos Kossuth University, Debrecen
(Bognár school; sugar-thiazolidine syntheses, refs 21-23). Not in `lab_ids.yml`; proposed id `debrecen_bognar`.
Linked to `gajda1990_extraction.md` only through the ligand synthesis method (Gajda's ref 3), not by authorship.

## 1. What this paper answers for the model

The case (M1c sec. 4.3) expected "C-2 inversion and C-4 epimerisation in TCAs" to give the P1 / P2 exchange speed
and the epimer share.

| model need | what Szilágyi 1979 gives | tag |
|---|---|---|
| C-2 epimer share at equilibrium (P2) | 2-aryl TCAs: ~1 : 1 (Table II: 45-55 % in Me2SO, 47-53 % pyridine, 46-54 % acetic acid; methyl ester 44-56 % in CDCl3 / C6D6); **no water, no sugar ring** | [M] |
| C-2 epimer exchange speed (P1 / P2) | **no rate.** Only: no NMR exchange broadening or coalescence of the two epimers' signals in Me2SO up to 140 °C, so the interconversion barrier is "much higher than 63 kJ/mol" on the NMR timescale | [M] / [Q] |
| C-2 inversion mechanism | during N-acetylation (Ac2O/water 100 °C, method A → 2,4-cis; Ac2O/pyridine 25 °C, method B → 2,4-trans) via ring-opened Schiff bases; no deuterium incorporation at C-2 in D2O/Ac2O at 100 °C, so not by C-2 deprotonation | [M] / [P] |
| C-4 epimerisation (cysteine racemisation) | only of the N-ACETYL acids in neat Ac2O at 90-100 °C (equilibrium trans : cis ≈ 4 : 5 within 10-20 min; full D at C-4); methyl esters unchanged ≥ 2 h at 90-100 °C, so it needs the free carboxyl (mixed anhydride, proposed) and does not open the ring | [M] / [P] |
| sugar-thiazolidine stereoselectivity | cited: sugar + cysteine condensations "sometimes remarkable stereoselectivity" (refs 21-22) | [C] |
| anything in water, any K or k | **none** | — |

## 2. Conditions

| item | value |
|---|---|
| Compounds | 1 2-phenyl-TCA (4R), 2 its 4S enantiomer, 3 2-p-tolyl-TCA (4R), 4 methyl 2-phenyl-TCA ester; N-acetyl 5-10, N-methoxycarbonyl 11, N-acetyl methyl esters 12-15; 16-17 N-benzylalanine references (Table I) |
| Syntheses | literature (refs 13, 42) for 1-3; 4 from (R)-cysteine methyl ester·HCl 10 mmol + benzaldehyde 10 mmol + pyridine in 16 mL methanol, overnight, room temperature |
| Acetylation, method A | 1 mmol TCA + 0.56 mL water + 0.50 mL Ac2O, steam bath (~100 °C) until homogeneous (5-10 min) |
| Acetylation, method B | 1 mmol TCA in 2.5 mL dry pyridine + 0.3 mL Ac2O, room temperature, 1-2 h |
| C-4 epimerisation | 1.0 g N-acetyl acid in 6 mL Ac2O, steam bath 20 min; NMR runs at 90 °C probe temperature |
| NMR | 1H 100 MHz (JEOL MH-100), 13C 25 MHz (Varian XL-100); solvents and temperatures per Table II (30-31 °C; 100 °C; up to 140 °C in Me2SO for 1 and 3) |
| Water / pH | **no aqueous kinetics; no pH** |

## 3. Tables re-typed

### Table II (p. 429), parent acids and their ester: epimer ratios (%) at 30-31 °C [M]

| compd | solvent | ratio (two sets) | H(2) ppm | H(4) ppm |
|---|---|---|---|---|
| 1 (2-phenyl-TCA) | (CD3)2SO | 55 / 45 | 5.70 / 5.49 | 4.25 / 3.92 |
| 1 | C5D5N | 47 / 53 | 6.25 / 5.80 | 4.58 / 4.32 |
| 1 | CH3COOH | 46 / 54 | 6.01 / 5.88 | 4.89 / 4.32 |
| 4 (methyl ester) | (CD3)2SO | 53 / 47 | 5.71 / 5.55 | 4.38 / 3.85 |
| 4 | C5D5N | 47 / 53 | 6.01 / 5.64 | 4.32 / 3.90 |
| 4 | CDCl3 | 44 / 56 | 5.77 / 5.52 | 4.16 / 3.95 |
| 4 | C6D6 | 47 / 53 | 5.73 / 5.43 | 3.90 / 3.75 |

Footnote a: "Epimer ratios in the case of 1 and 4, and conformer population ratios in the other cases." The N-acetyl
compounds (5, 8, 9, 10, 12-14: 54-91 % major set) are amide-rotamer populations, which coalesce below 100 °C [M/Q];
not transcribed further.

### Table III (p. 430) and Experimental (p. 431): C-4 epimerisation in Ac2O [M]

- 5 (2R,4R) → 5 ([α]D +129°) and 7 (2R,4S, +304°); 8 (2S,4R) → 8 (−306°) and 6 (2S,4S, −121°).
- 5 or 8 in Ac2O at 90 °C, 20 min: two H(4) signals in **4 : 5** intensity (trans 8 : cis 5); the same mixture from
  either side.
- Text: "approximately 1:1" mixture after "a short time (10 min)" at 100 °C.
- 12 and 14 (methyl esters) unchanged after ≥ 2 h at 90-100 °C in Ac2O.
- 1-d2 in Ac2O-d6, 90 °C: 75 % D at C-4, ≈ complete given the 80 % isotopic purity of the start.

## 4. Arithmetic (mine) [Z]

**1. The NMR bound is not a pot-timescale bound.** "No coalescence up to 140 °C" means the epimer interconversion of
1 and 3 in Me2SO is slower than ~10¹-10² s⁻¹ (signal separations of ~0.2 ppm at 100 MHz ≈ 20 Hz). That is
an upper bound of order 10³-10⁴ min⁻¹, far above any pot rate. It says nothing about whether exchange is minutes or
days in water.

**2. Epimer share.** 44-56 % in every solvent: for these 2-aryl rings the 2,4-cis / 2,4-trans free-energy difference
is ≤ 0.6 kJ/mol at 30 °C (RT ln(56/44)). TTCA's 2R : 2S = 4.06 at 37 °C (Radomski 1989) is 3.6 kJ/mol: the
polyhydroxyalkyl C-2 is much more selective than an aryl one.

## 5. What it can and cannot test in the model

- **Can:** show that C-2 epimers of cysteine thiazolidines form non-stereospecifically and survive heating to 140 °C
  in an aprotic solvent without visible exchange broadening; that C-4 (the cysteine stereocentre) epimerises only
  through an activated carboxyl, so cysteine racemisation via the ring is not expected in water pots [P, transfer].
- **Cannot:** give any rate or equilibrium constant in water; any sugar-ring epimer ratio; any pH dependence; any
  P1a open fraction.

## 6. Conflicts with existing records

| record | what it says | what Szilágyi says | verdict |
|---|---|---|---|
| `radomski1989_extraction.md` (TTCA 2R : 2S = 4.06, D2O 37 °C; relaxation in minutes) | sugar ring selective, fast exchange in water | aryl rings ~1 : 1, no exchange broadening in Me2SO | **no conflict**: different class and solvent; the NMR bound is consistent with minutes-scale exchange |
| `nagasawa1981_extraction.md` (C-2 epimerisation hours in neutral water, < 1 h in 1 N HCl) | acid-catalysed ring-chain step | C-2 inversion proposed via Schiff base during acylation | consistent in mechanism; no rate here |
| `pesek1975_extraction.md` (aliphatic TCAs stable for days in neutral D2O) | stable rings | 2-aryl TCAs show no NMR change to 140 °C in Me2SO | consistent (stability, not exchange) |
| Parthasarathy 1976 (this paper's ref 24; not on disk, cited by no dossier) | p-tolualdehyde + cysteine gives one epimer, two "very stable conformers" | they are two C-2 epimers, ~1 : 1 | recorded only so that no later reading takes ref 24's single-epimer claim at face value |
| `gajda1990_extraction.md` (sugar ligands 1 : 1 after recrystallisation) | 1 : 1 | 1 : 1 for aryl rings | consistent for crystals; neither is a water equilibrium for TTCA |

## 7. Flags

1. **[!] No kinetics in water.** The queue title suggested "ring-chain rates"; the paper reports none.
2. All exchange experiments are in acetic anhydride (± water, ± pyridine) or deuterated organic solvents.
3. Mechanisms (Schemes I-II) are proposals [P], supported by D-incorporation and product analysis.
4. The text layer drops Table II's solvent column; values above are from the render.

## M1c role (frozen bd641ae8)

**Role fixed by the case (sec. 4.3):** band source for the P1 / P2 exchange speed and the epimer share.

**Band it can supply:**
- **P2 (epimer share):** nothing for TTCA. The 44-56 % shares are for 2-ARYL rings in organic solvents at 30 °C;
  the band 0-0.20 of TTCA (Radomski 1989, water, 37 °C, 2R : 2S = 4.06) is unchanged.
- **P1a / P2 exchange speed:** no rate. The only number is an NMR-timescale bound (no coalescence to 140 °C in
  Me2SO), which is weaker by orders of magnitude than the case's existing ≥ 2 × 10⁻³ min⁻¹ class bound (Nagasawa
  1981) and does not narrow it.
- **A qualitative point the case may use:** C-4 (cysteine) epimerisation needs an activated carboxyl and does not
  open the ring; cysteine racemisation is not a hidden P-member in water pots.
- **Transfer:** distant (aryl C-2, organic solvents). **Not decisive**, and it measures nothing the case would call
  decisive.

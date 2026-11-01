# Doorn & Petersen 2002 — EXTRACTION (4-hydroxy-2-nonenal (4HNE) and 4-oxo-2-nonenal (4ONE) with N-acetyl amino-acid derivatives and peptides: second-order rate constants, 50 mM phosphate pH 7.4, ~23 °C, by UV loss of the enal; peptide adducts by MALDI-TOF-MS / PSD, 37 °C. Thiols are N-acetylcysteamine, GSH and a Cys peptide; NO free cysteine)

### SELECTIVITY NUMBERS, ONE pH, ONE TEMPERATURE. Thiol ≫ His > Lys (> Arg for 4ONE): for 4HNE, thiols 1.2-6.8 M⁻¹ s⁻¹ vs His 2.1 × 10⁻³ and Lys 1.3 × 10⁻³; 4ONE reacts 110-150x faster than 4HNE with thiols (145-887 M⁻¹ s⁻¹). Forward rates only: no reverse rate, no K, no pH law, no Ea.

**Read 2026-11-01 from `data/articles/doorn2002.pdf`** (6 pp., Chem. Res. Toxicol. 15, 1445-1450; born-digital, clean
text layer). Prose and Table 2 from `pdftotext -layout`; **Table 2 and Figs. 5-6 checked against 130 dpi renders**,
pp. 1446 against an 80 dpi render (scratchpad `dos6/doo-4.png`, `doo1-2.png`). Figure reads ±5 %.

**Legend.** **[V]** verbatim · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative claim ·
**[C]** cited · **[P]** proposed · **[Z]** derived by me · **[!]** flag.

## 1. What this paper answers for the model

| model need | what Doorn 2002 gives | tag |
|---|---|---|
| forward k, thiol + enal | 4HNE: N-acetylcysteamine 1.21 ± 0.03, GSH 1.33 ± 0.05, Ac-RRWWCR-amide 6.79 ± 0.65 M⁻¹ s⁻¹; 4ONE: 186 ± 29, 145 ± 10, 887 ± 88 M⁻¹ s⁻¹ (Table 2, mean ± SE, 23 °C, pH 7.4) | [M] |
| forward k, amine / imidazole + enal | 4HNE: Ac-Lys-amide 1.33 ± 0.08 × 10⁻³, Ac-histamine 2.14 ± 0.31 × 10⁻³; 4ONE: 7.46 ± 0.61 × 10⁻³, 2.21 ± 0.17 × 10⁻², Ac-Arg-amide 6.11 ± 0.85 × 10⁻⁴ M⁻¹ s⁻¹ | [M] |
| reverse k, K | **none** | — |
| pH law, buffer catalysis | **none** (one pH, one buffer); the Arg-rich Cys peptide is 5-6x faster, ascribed to a lowered thiol pKa | [M] / [P] |
| Ea | **none** (~23 °C kinetics; 37 °C adduction for MS only) | — |
| stoichiometry | 1:1 Michael adducts (+156 Da 4HNE, +154 Da 4ONE); 4ONE + GSH is a 1:1 bimolecular titration (Fig. 6); bis-adducts on peptides with Arg + Cys / Lys | [M] |

**Answers.**

1. **MFT loss law: no bound** (confidence ~95 %). These are 23 °C forward constants for lipid enals with alkyl thiols at
   pH 7.4. What transfers [Z, assumption: MFT reacts as an alkyl thiolate of similar pKa]: (i) thiol capture by a lipid
   enal is second order in the electrophile (A2's form); (ii) the electrophile's structure moves the constant by 10²
   (4-oxo vs 4-hydroxy), the same scale as Esterbauer's structure series; (iii) λ = k[E]: at 23 °C, 1 µM 4ONE gives
   λ ≈ 0.009 min⁻¹ for any thiol at pH 7.4, 1 µM 4HNE ≈ 8 × 10⁻⁵ min⁻¹.
2. **Selectivity** [M → Z]: per mole, thiols beat Lys amines ~900x (4HNE) and ~2.5 × 10⁴x (4ONE). In a lipid-bearing
   pot with protein, free cysteine and MFT compete for unsaturated lipid aldehydes mainly with other THIOLS, not with
   amines: cysteine at tens of mM shields MFT (the A3-type competition of `esterbauer1975_extraction.md`), protein
   lysine does not.
3. **CF-031 (lipid lowers vs raises meaty thiols)**: supports, in kind, "lipid aldehydes as thiol sinks" for 4HNE / 4ONE
   type oxidation products; does not cover hexanal or other saturated aldehydes (not Michael acceptors). LEAVES the
   conflict open (no pot data).
4. **`k_ttca_cys`, CF-076 E1 / E2 / E4**: nothing. **E3 (second Cys sink)**: nothing directly (no free cysteine, no
   dicarbonyl).

## 2. Identity

| field | value |
|---|---|
| Title | "Covalent Modification of Amino Acid Nucleophiles by the Lipid Peroxidation Products 4-Hydroxy-2-nonenal and 4-Oxo-2-nonenal" |
| Authors | Jonathan A. Doorn, Dennis R. Petersen* — School of Pharmacy, University of Colorado Health Sciences Center, Denver |
| Venue | *Chem. Res. Toxicol.* **2002**, 15(11), 1445-1450; received August 1, 2002 |
| DOI | **10.1021/tx025590o** (printed p. 1445; Crossref 2026-11-01: same title, Doorn & Petersen, 15(11) 1445-1450). Matches the brief |
| Funding | NIH / NIEHS / NIAAA grants (acknowledgement) |

## 3. Conditions

- **Kinetics** [M]: 50 mM sodium phosphate pH 7.4, ~23 °C; nucleophile ≥ 10x over 4HNE / 4ONE (pseudo-first order);
  enal loss at 224 nm (4HNE) / 228 nm (4ONE), plate reader or diode array; thiol alternative: 4HNE > 10x over thiol,
  quench with excess DTNB, 412 nm. k′ from ln(% enal) vs t; k from the slope of k′ vs [nucleophile]; ≥ 3 independent
  experiments per k′. NAC here = **N-acetylcysteamine** (not N-acetylcysteine) [!].
- **MS adduction** [M]: peptides 10 µM + 200 µM 4HNE / 4ONE, 10 mM tricine pH 7.4, 37 °C, 16-18 h; MALDI-TOF-MS and PSD.
- 4HNE from its dimethyl acetal; 4ONE by Dess-Martin oxidation of 4HNE (UV λmax 228 nm, NMR, TLC).
- O₂ / headspace: not stated.

## 4. Quantitative results

**Table 2** (23 °C, pH 7.4; k in M⁻¹ s⁻¹, mean ± SE) [M]:

| nucleophile | 4ONE | 4HNE | k4ONE / k4HNE |
|---|---|---|---|
| Ac-Lys-amide | 7.46 (± 0.610) × 10⁻³ | 1.33 (± 0.0830) × 10⁻³ | 5.61 ± 0.577 |
| Ac-histamine | 2.21 (± 0.170) × 10⁻² | 2.14 (± 0.312) × 10⁻³ | 10.3 ± 1.70 |
| Ac-Arg-amide | 6.11 (± 0.850) × 10⁻⁴ | ND (no adduct by MS) | ND |
| Ac-cysteamine | 186 ± 29.3 | 1.21 ± 0.0268 | 153 ± 24.3 |
| GSH | 145 ± 10.1 | 1.33 ± 0.0500 | 109 ± 8.63 |
| Ac-RRWWCR-amide | 887 ± 88.3 | 6.79 ± 0.646 | 130 ± 18.0 |

Footnote d: GSH + 4HNE "in agreement" with Esterbauer's 1.09 M⁻¹ s⁻¹ [C].

**Fig. 5** [F]: 4HNE + N-acetylcysteamine 0.5-5 mM, k′ 1.1 × 10⁻³ → 6.4 × 10⁻³ s⁻¹ (slope ≈ 1.2 M⁻¹ s⁻¹); 4ONE,
0.1-0.5 mM, k′ 0.022 → 0.102 s⁻¹ (slope ≈ 195). **Fig. 6** [F]: GSH and 4ONE 100 µM each, ~45 µM left at 60 s, ~37 at
80 s (second-order t½ ≈ 60 s → k ≈ 170 M⁻¹ s⁻¹ [Z], consistent with Table 2). Text: with 5 mM GSH, 4ONE t½ < 1 s [M].

## 5. Mechanism evidence and its kind

- Colorado laboratory (Petersen). Kind 6 within one study (nucleophile identity, electrophile identity at fixed pH and
  T); structure by MS mass increments and PSD fragmentation (Michael adducts; 4ONE adduct predicted to keep two
  carbonyls) [M/P]. Reversibility not tested.

## 6. Bearing on open questions

| item | verdict | why |
|---|---|---|
| MFT loss law, A2 (electrophile-pool second order) | SUPPORTS the form in kind | second order in enal, structure-dependent 10² |
| thiol-rivals A3 (cysteine shields) | SUPPORTS competition (thiol vs thiol), says amines do not shield | thiol / Lys ~10³-10⁴ |
| CF-031 | SUPPORTS feasibility of unsaturated lipid aldehydes as thiol sinks; leaves the conflict open | no pot data |
| CF-076, `k_ttca_cys` | LEAVES UNCHANGED | |

## 7. Disagreements found

1. With `esterbauer1975_extraction.md`: GSH + 4HNE 1.33 (23 °C, 50 mM phosphate) vs 1.09 (20 °C, 66 mM phosphate):
   agreement within 22 %, a passed cross-lab check, not a disagreement.
2. Text vs table: the Discussion says "pH 7.4, 37 °C" for the adduct-stability statement; the kinetics are at 23 °C
   (Table 2 footnote). Not an inconsistency once the two experiments are separated.

## 8. Ledger-ready measurements

1. **Table 2** (all six rows × two enals): printed, SE given; the only thiol-vs-amine selectivity constants on disk for
   lipid enals.

# Davis, Layland, Page, Martin & More O'Ferrall 1991 (Part 2) — EXTRACTION
DOI: 10.1039/p29910001225 (not printed in the PDF; Crossref 2026-11-05: "Thiazolidine ring opening in penicillin derivatives. Part 2. Enamine formation", Davis, Layland, Page, Martin, O'Ferrall, J. Chem. Soc. Perkin Trans. 2 1991, 1225. **Matches the PDF; the brief's author list "Davis, Jones & Page" is Part 1's** (Flag 1).)
### HYDROXIDE-CATALYSED THIAZOLIDINE RING OPENING BY ELIMINATION TO AN ENAMINE (across C-6/C-5 = across the exocyclic C-1' / ring C-2 in thiazolidine terms), first order in OH⁻, for the methyl ester and two amides of benzylpenicilloate and (slowly) for the penicilloate dianion, 30 °C, I 1.0 M KCl. Ester: k1 = 8.09 × 10⁻² M⁻¹ s⁻¹ (opening), k−1 = 6.42 × 10⁻² s⁻¹ (pH-independent re-closure), K = 1.26 M⁻¹; amides 0.115 / 0.15; dianion 4.50 × 10⁻⁴ M⁻¹ s⁻¹. Solvent isotope effect 7.5 on re-closure. **The structural analogue of a thiazolidine → Amadori-type (enaminol) exit, with the C-H made acidic by an adjacent ester or amide; no constant transfers to TTCA.**

**Read 2026-11-05 from `data/articles/davis1991b.pdf`** (5 pp., Perkin Trans. 2 pp. 1225-1229; RSC scan with a text
layer). Prose from `pdftotext -layout`; **the rate constants of p. 1227 and the dianion constant of p. 1228 re-read
from 200 dpi renders** (scratchpad `img/davb3r.png`, `davb4l.png`), because the text layer garbles their exponents.
Two exponents on p. 1228 ("ca. 0.01 for a 10⁻? mol dm⁻³ solution … in 10⁻? mol dm⁻³ sodium hydroxide") are lost in
the render as well and are not used.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Huddersfield Polytechnic (Page) with University College Dublin (Martin, More O'Ferrall). One lab id
with `davis1991_extraction.md` and, through More O'Ferrall, linked to `canle1996_extraction.md` (proposed ids
`huddersfield_page`, `ucd_more_oferrall`; merge rule 2 would join them).

## 1. What this paper answers for the model

The case lists Part 2 with Part 1 as a band source for "the pH shape of ring opening (R-rel)". Part 2 measures a
DIFFERENT opening: not the unimolecular C-S cleavage to an iminium (Part 1), but a base-catalysed elimination that
opens the ring directly to an enamine.

| model need | what Davis 1991 Part 2 gives | tag |
|---|---|---|
| a thiazolidine exit that is first order in OH⁻ | enamine formation, rate ∝ [OH⁻] (0.02-0.3 M NaOH; and "below pH 11" for the ester, same second-order constant) | [M] |
| its rate constants, 30 °C | ester 4: **k1 = 8.09 × 10⁻² M⁻¹ s⁻¹**; amide 7: **0.115**; N-ethylamide 8: **0.15**; penicilloate dianion 6: **4.50 × 10⁻⁴** (only > 1 M NaOH) | [M] |
| re-closure and equilibrium | k−1 = **6.42 × 10⁻² s⁻¹** (pH-independent); K = k1/k−1 = **1.26 M⁻¹**; enamine ≤ ~10 % of the ester under the run conditions | [M] |
| competing step | methyl-ester hydrolysis k2 = 0.327 M⁻¹ s⁻¹ | [M] |
| mechanism markers | k−1(H2O)/k−1(D2O) = 7.5; inverse 0.88 on k1; 0.97 on k2; concerted E2-type elimination preferred [P] | [M] / [P] |
| neutral-pH opening | "slow reversible unimolecular" opening to an iminium, faster than ester hydrolysis over pH 4-11 | [Q] |
| anything for a sugar ring, any T ≠ 30 °C | **none** | — |

## 2. Conditions

| item | value |
|---|---|
| Substrates | (5R,6R)-methyl benzylpenicilloate 4 (made as in Part 1 / ref 7); α-amide 7 and α-N-ethylamide 8 of benzylpenicilloic acid (from benzylpenicillin + NH3 / EtNH2; analyses, NMR, IR, MS printed); penicilloic acid 6 |
| Medium | freshly boiled glass-distilled water, **I = 1.0 M KCl**; NaOH 0.02-0.10 M (Fig. 1) and to 0.3 M (Figs. 2-3); buffers for pH < 11 not named |
| Temperature | **30 ± 0.1 °C** (cell block; stopped-flow bath 30 ± 0.1 °C) |
| pH | Philips PW 9409 + Russell CEL electrode calibrated at 30 °C; checked before and after each run, drift ≤ 0.03 |
| Detection | UV 280-290 nm (Gilford 2600; Nortech SF3A stopped flow for t½ < 2 s); consecutive first-order fits (λ1, λ2); NMR in NaOD/D2O |
| Substrate concentration | 1.0 × 10⁻? M (Fig. 1 caption; exponent lost) [!] |

## 3. Numbers re-typed (pp. 1227-1228) [M]

| quantity | value | substrate |
|---|---|---|
| k−1 (enamine → ring, pH-independent) | 6.42 × 10⁻² s⁻¹ | ester 4 |
| k1 (OH⁻-catalysed ring opening to enamine) | 8.09 × 10⁻² dm³ mol⁻¹ s⁻¹ | ester 4 |
| k2 (ester hydrolysis) | 0.327 dm³ mol⁻¹ s⁻¹ | ester 4 |
| ε285 of the penamaldate | 17 200 dm³ mol⁻¹ cm⁻¹ | ester 4 |
| K = k1/k−1 | 1.26 dm³ mol⁻¹ | ester 4 |
| k (chromophore, > 1 M NaOH) | 4.50 × 10⁻⁴ dm³ mol⁻¹ s⁻¹ ("ca. 200 times slower" than the ester) | penicilloate 6 |
| k (chromophore, 290 nm, corrected for amide ionisation) | 0.115 / 0.15 dm³ mol⁻¹ s⁻¹ | amides 7 / 8 |
| solvent KIE | k−1: 7.5; k1: 0.88; k2: 0.97; λ1λ2 slope in D2O 5.4x lower | ester 4 |
| thiol pKa of the open enamine-thiol 9 (assumed) | 11.0 → equilibrium 9 : 4 ≈ 10⁻³ | [P] |
| comparison, Part 1 | C-5 epimerisation of 6 by OH⁻: 7.6 × 10⁻⁴ dm³ mol⁻¹ s⁻¹ | [C, Part 1] |

No errors printed. Absorbance change independent of pH below pH 10.2, decreasing above [Q].

## 4. Arithmetic (mine) [Z]

**1. What the elimination is, in TTCA terms.** In penicilloate the ring opens by removing H-6 (on the exocyclic
carbon next to ring C-5) while C-5-S breaks, giving C-5=C-6 with N on C-5: an enamine. In TTCA the same bond pattern
is H on C-1' (the sugar's C-2, bearing OH) next to ring C-2 (S, N): the analogous elimination gives a 1,2-enaminol,
i.e. the Amadori route (TTCA → CAMA) without a free Schiff base. **This is a mechanistic precedent for a
hydroxide-catalysed (n = 1) TTCA → CAMA exit (the case's X-cama).** It is not a rate for it: the C-H here is
acidified by an α-ester or α-amide (k1 0.08-0.15 M⁻¹ s⁻¹); with an α-carboxylate it is 200x slower; TTCA's C-1'-H
has only an α-hydroxyl.

**2. Order of magnitude at the pots' pH (heavy transfer, for scale only).** At 30 °C and pH 8 ([OH⁻] ≈ 1.5 × 10⁻⁶ M),
even the ester's k1 gives 1.2 × 10⁻⁷ s⁻¹ = 7 × 10⁻⁶ min⁻¹; the dianion's 4 × 10⁻⁸ min⁻¹. The fed-TTCA loss the case's
part B describes is 0.0164 min⁻¹ at label pH 8, 100 °C (Zhai 2021 Fig. 4b). With pKw ≈ 12.3 at 100 °C and in-situ pH
8, [OH⁻] ≈ 5 × 10⁻⁵ M, so an OH⁻ route would need k_OH ≈ 5 M⁻¹ s⁻¹ at 100 °C: ~60x the activated ester's 30 °C value
(an Ea near 50 kJ/mol would give ~30-60x from 30 to 100 °C) and ~10⁴x the dianion's. For TTCA's unactivated C-H that
looks too fast; this is arithmetic on a transfer, not a test.

**3. Equilibrium position.** K = 1.26 M⁻¹: at 0.1 M OH⁻ about 11 % of the ester sits as enamine at equilibrium
(K[OH⁻]/(1 + K[OH⁻])), matching "never exceeds ca. 10 %".

## 5. What it can and cannot test in the model

- **Can:** show that a thiazolidine can open by a base-catalysed elimination that is exactly first order in OH⁻ and
  reversible (re-closure ~0.06 s⁻¹ at 30 °C), and that its rate scales with the acidity of the exocyclic C-H (ester ≈
  amide ≫ carboxylate). A form precedent for a one-hydroxide (n = 1) exit to an Amadori-type product.
- **Cannot:** give a magnitude for TTCA (activated vs unactivated C-H; 30 °C; ≥ 0.02 M NaOH for the measured rows);
  tell X-cama from X-deg; give any release (hydrolysis) rate.

## 6. Conflicts with existing records

| record | what it says | what Part 2 says | verdict |
|---|---|---|---|
| `davis1991_extraction.md` (Part 1: ring opening flat pH 6-12) | no OH⁻ term below pH 12.5 for the dianion | dianion's enamine route 4.50 × 10⁻⁴ M⁻¹ s⁻¹, only > 1 M NaOH | consistent (same lab; the OH⁻ routes need an acidic C-H to matter) |
| `pesek1975_extraction.md` (aliphatic TCAs decompose only in > 3 M NaOH; OH⁻-catalysed ring-chain k1 25 M⁻¹ s⁻¹ for aryl TCAs) | base opens at extremes | base opening slow unless C-H activated | consistent |
| de Roos 2005 (in the case: "At pH < 7, the route to Amadori compounds is preferred") and Zhai 2020 (TTCA → CAMA favoured at pH 7.5-9.5, low water) | opposite qualitative directions | an OH⁻-catalysed elimination route exists in principle | no data to choose; the paper supplies a mechanism, not a direction for TTCA |

## 7. Flags

1. **[!] Authors.** The brief lists "Davis, Jones & Page" for both parts; Part 2's byline (and Crossref) is Davis,
   Layland, Page, Martin & More O'Ferrall. Same volume, consecutive pages 1225-1229.
2. **[!] Measured rows are mostly in 0.02-0.3 M NaOH**, pH ≈ 12.3-13.5; the "below pH 11" statements give no table.
3. The penamaldate ε and K depend on the λ1/λ2 decomposition (Scheme 3), which the authors fit three ways with
   "very good agreement" [Q]; no errors.
4. 30 °C, I 1.0 M only.

## M1c role (frozen bd641ae8)

**Role fixed by the case (sec. 4.3):** band source for the pH shape of ring opening (R-rel), jointly with Part 1.

**Band it can supply:**
- **For R-rel's n_r:** nothing new on release; Part 2 does not measure hydrolysis. Part 1 carries that band
  (n_r ≤ 0 over pH 4-8).
- **For part B's slope n (the exit), as a FORM precedent only:** a measured thiazolidine ring opening to an
  Amadori-type enamine that is exactly first order in OH⁻ (n = 1 dex per pH unit), i.e. the upper end of the case's
  0-1 sensitivity, and the "one hydroxide" form the case names. Magnitudes (30 °C, I 1.0 M): 4.50 × 10⁻⁴ (α-carboxylate)
  to 0.15 M⁻¹ s⁻¹ (α-amide); re-closure 6.42 × 10⁻² s⁻¹. The case's n band (0.30-0.74, Martins & van Boekel) is not
  changed; the paper does not say TTCA follows this route.
- **Transfer:** distant for magnitude (activated C-H, 30 °C, strong base), same-class for form.
- **Not decisive.** It measures nothing the case would call decisive.

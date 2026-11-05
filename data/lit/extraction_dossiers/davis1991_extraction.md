# Davis, Jones & Page 1991 (Part 1) — EXTRACTION
DOI: 10.1039/p29910001219 (not printed in the PDF; Crossref 2026-11-05: "Thiazolidine ring opening in penicillin derivatives. Part 1. Imine formation", Davis, Jones, Page, J. Chem. Soc. Perkin Trans. 2 1991, 1219. **Matches the brief and the PDF.** Part 2 is `davis1991b_extraction.md`.)
### A FULL pH-RATE PROFILE (pH 1-13.5, 30 °C, I 1.0 M KCl) FOR THIAZOLIDINE RING OPENING in an aliphatic, penicillamine-type ring, read through C-5 (= thiazolidine C-2) epimerisation of benzylpenicilloate and its mono- and di-methyl esters. Ring opening is **pH-independent from pH 6 to 12** (k_obs 5.84 × 10⁻⁵ s⁻¹, unimolecular C-S cleavage ≥ 4.96 × 10⁻⁵ s⁻¹), **rises 6.8x on the acid side** of the ring-N pKa (5.14 / 4.59 for the two epimers) and becomes **hydroxide-dependent only above pH 12.5** (7.60 × 10⁻⁴ M⁻¹ s⁻¹). Ring closure > 10⁸ s⁻¹. Hydrolytic release (to penicillamine + aldehyde) is ≥ 20x slower than opening below pH 7 and not seen above pH 8.

**Read 2026-11-05 from `data/articles/davis1991.pdf`** (5 pp., Perkin Trans. 2 pp. 1219-1223; RSC scan with a text
layer). Prose from `pdftotext -layout`; **Table 1 (p. 1223), Fig. 1 and eqn. (3) (p. 1221) read from 220 dpi renders**
(scratchpad `img/dav-3.png`, `dav5tab.png`), because the text layer garbles every exponent in Table 1. Fig. 1 read
error ± 0.05 log unit, ± 0.2 pH [F].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Chemical and Physical Sciences, Huddersfield Polytechnic (Page; physical organic,
β-lactams). Not in `lab_ids.yml`; proposed id `huddersfield_page`. Part 2 adds More O'Ferrall (Dublin), linking this
lab to `canle1996_extraction.md` under merge rule 2.

## 1. What this paper answers for the model

The case (M1c sec. 4.3) wants "the pH shape of ring opening (R-rel)": is the ring-chain step, or the release it
gates, base-favoured?

| model need | what Davis 1991 Part 1 gives | tag |
|---|---|---|
| ring-opening (ring ⇌ iminium-thiolate) rate vs pH, aliphatic C-2 | pH-rate profile pH 1-13.5 (Fig. 1); plateau pH 6-12; acid limb; OH⁻ limb above 12.5 | [M] / [F] |
| unimolecular C-S cleavage, neutral ring | **k_f⁰ = 4.96 × 10⁻⁵ s⁻¹** (30 °C; a lower limit on ring opening) | [M] |
| same, N-protonated ring (kinetically S-protonated) | k_f^H = 2.45 × 10⁻⁴ s⁻¹ | [M] |
| ring-N pKa (penicilloate) | 5.14 (5R), 4.59 (5S); monoester 3.80 / 3.35; diester 1.18 / 1.51 | [M] |
| epimer equilibrium at C-2 (thiazolidine numbering) | K⁰ = 5.67 (5S : 5R = 85 : 15, neutral); K^H = 1.60 (60 : 40, protonated) | [M] |
| ring closure | > 10⁸ s⁻¹ (bound, from pH-independence at pH 11 and estimated proton-transfer rates) | [Z by the authors] |
| hydroxide-catalysed opening (to the neutral imine) | 7.60 × 10⁻⁴ M⁻¹ s⁻¹ (pH > 12.5) | [M] |
| hydrolytic release (aminothiol + aldehyde) | not measured; "at least 20 times" slower than epimerisation below pH 7; no degradation products above pH 8 "after several days" | [M] / [Q] |
| buffer catalysis | not reported (buffers not named) | [!] |

## 2. Conditions

| item | value |
|---|---|
| Substrates | (3S,5R,6R)-benzylpenicilloate 5 (from benzylpenicillin by B. cereus β-lactamase; epimerically pure); its α-monomethyl ester 8 and dimethyl ester 9. Ring = 5,5-dimethyl (penicillamine) thiazolidine-4-carboxylate whose C-2 carries CH(NHCOCH2Ph)CO2⁻ / CO2Me |
| Solvent / ionic strength | freshly boiled glass-distilled water; **I = 1.0 M KCl** "except where otherwise indicated"; buffers from AnalaR chemicals, **identities and concentrations not stated** |
| Temperature | **30 °C** |
| pH | Philips PW 9409, Russell CEL combination electrode calibrated against standard buffers at 30 °C (activity scale) |
| Analysis | HPLC (RP C18 or PS-DVB, 254 nm), 270 MHz 1H NMR (D2O), polarimetry; first-order fits by non-linear least squares; very slow runs by initial slopes |
| Concentrations | not printed for the kinetic runs |

## 3. Tables and figures re-typed

### Table 1 (p. 1223). C-5 epimerisation, water, 30 °C [M]

| derivative | k_obs⁰ / s⁻¹ | k_f⁰ / s⁻¹ | k_r⁰ / s⁻¹ | K⁰ | k_obs^H / s⁻¹ | k_f^H / s⁻¹ | k_r^H / s⁻¹ | K^H | pKa^R | pKa^S |
|---|---|---|---|---|---|---|---|---|---|---|
| diacid 5 | 5.84 × 10⁻⁵ | 4.96 × 10⁻⁵ | 8.74 × 10⁻⁶ | 5.67 | 3.98 × 10⁻⁴ | 2.45 × 10⁻⁴ | 1.53 × 10⁻⁴ | 1.60 | 5.14 | 4.59 |
| monoester 8 | 2.82 × 10⁻⁶ | 1.12 × 10⁻⁶ | 1.70 × 10⁻⁶ | 0.66 | 2.37 × 10⁻⁴ | 1.42 × 10⁻⁴ | 9.48 × 10⁻⁵ | 1.50 | 3.80 | 3.35 |
| diester 9 | 3.53 × 10⁻⁸ | 1.39 × 10⁻⁸ | 2.14 × 10⁻⁸ | 0.65 | 3.16 × 10⁻⁴ | 1.84 × 10⁻⁴ | 1.31 × 10⁻⁴ | 1.40 | 1.18 | 1.51 |

No errors printed. ⁰ = unprotonated thiazolidine, H = N-protonated; f = 5R → 5S. Eqn. (3): k_f⁰ = k_obs⁰ / (1 + 1/K⁰);
eqn. (4) combines both forms with the fractions α_R, α_S of unprotonated epimers.

### Text constants (pp. 1220-1222) [M]

- pH-independent epimerisation of 5: 5.84 × 10⁻⁵ s⁻¹; no D incorporation at 5-H or 6-H in D2O; < 2 % 6S epimers in
  10 h; no enamine by UV.
- "minimum rate constant for unimolecular cleavage of the C-S bond" of 5 at 30 °C: 4.96 × 10⁻⁵ s⁻¹.
- Diester 9: 1.7 × 10³ times slower than 5 at the plateau; monoester 8: 21 times slower.
- Thiol pKa of the open iminium-thiol estimated 7.9, iminium pKa 5.3 (ref 20); ring closure "must be greater than
  10⁸ s⁻¹" [authors' estimate].
- Base-catalysed C-5 epimerisation of 5 above pH 12.5: second order, **7.60 × 10⁻⁴ dm³ mol⁻¹ s⁻¹**; proposed
  concerted E2-type opening to the neutral imine [P].
- Enamine 7 chromophore (280 nm) only above 1 M NaOH: 4.50 × 10⁻⁴ dm³ mol⁻¹ s⁻¹.
- C-6 epimerisation of benzylpenicillin (with D exchange) by OH⁻: ca. 6 × 10⁻³ dm³ mol⁻¹ s⁻¹ (not accurately
  determined).
- Diester 9 at pH 1 is 10⁴ times faster than at pH 7.
- Degradation below pH 7 (to benzylpenilloic acid, then D-penicillamine + benzylpenilloaldehyde): epimerisation
  "at least 20 times greater" than degradation (< 5 % products during epimerisation). Above pH 8 penicilloate "very
  stable", no degradation products by HPLC after several days.

### Fig. 1 (p. 1221), diacid 5, log k_obs (s⁻¹) vs pH [F]

≈ −3.40 at pH 1-2.5; −3.85 at pH 5; −4.15 to −4.25 from pH 6.3 to 11.3; −4.0 at pH 12.8; −3.6 at 13.2; ≈ −3.3 at 13.4.
The line is eqn. (4) with Table 1.

## 4. Arithmetic (mine) [Z]

**1. In model units.** Plateau k_obs⁰ = 3.5 × 10⁻³ min⁻¹, k_f⁰ = 3.0 × 10⁻³ min⁻¹; protonated k_obs^H = 2.4 × 10⁻²
min⁻¹ (30 °C). Acid limb / plateau = 6.8x = +0.83 dex across the ring-N pKa (≈ 4.6-5.1).

**2. Where hydroxide matters.** k_OH [OH⁻] = k_obs⁰ at [OH⁻] = 5.84 × 10⁻⁵ / 7.60 × 10⁻⁴ = 0.077 M, pH ≈ 12.7-12.9 at
30 °C. At pH 8 the hydroxide term is ~10⁻⁵ of the plateau.

**3. Release bound.** Below pH 7, degradation ≤ k_epimerisation / 20 ≈ ≤ 1.8 × 10⁻⁴ (plateau) to ≤ 1.2 × 10⁻³ min⁻¹
(acid limb) at 30 °C. Above pH 8, "several days" with no product: of order < 10⁻⁵ min⁻¹ if the HPLC would see ~2 %
in ~3 days (my assumed detection; the paper gives none). Release is thus faster at pH < 7 than at pH > 8.

## 5. What it can and cannot test in the model

- **Can:** give the class pH-SHAPE of thiazolidine ring opening for an aliphatic C-2 carboxylate ring: flat over pH
  6-12, faster on the acid side of the ring-N pKa, hydroxide-dependent only at pH > 12.5; and the class statement that
  opening (exchange) is much faster than hydrolytic release, which is acid-favoured too.
- **Cannot:** give a constant for TTCA (C-2 carries an acylamino-carboxylate, not a polyol; gem-dimethyl at C-5);
  any temperature dependence (30 °C only); any buffer term; the release rate itself.

## 6. Conflicts with existing records

| record | what it says | what Davis 1991 says | verdict |
|---|---|---|---|
| `nagasawa1981_extraction.md` (penicillamine + acetaldehyde ring: C-2 epimerisation ~24 h neutral, < 1 h in 1 N HCl, stops in Na2CO3) | acid-catalysed, base-slowed | acid limb 6.8x; flat to pH 12 | **consistent in direction**; Davis's plateau (t½ ~3 h at 30 °C) is faster than Nagasawa's ~24 h, and Davis sees no slowing in base up to pH 12 |
| `fife1991_extraction.md` (ring opening fast, iminium hydrolysis rate-determining and general-base catalysed) | release is base/buffer-catalysed for aryl rings | release here is acid-side; none above pH 8 | not a direct conflict (different C-2 class and pH range), but **the aliphatic class does not show Fife's base-favoured release** |
| `karpinski1990_extraction.md` (AThz release 1.8 × 10⁻⁴ → 4.8 × 10⁻⁵ min⁻¹ from pH 6.0 to 7.4, 25 °C) | release falls with pH | release ≥ 20x below opening at pH < 7, absent at pH > 8 | consistent |
| `pesek1975_extraction.md` (aliphatic TCAs stable in neutral D2O, decompose only in > 3 M NaOH) | base opens only at extremes | OH⁻ term only above pH 12.5 | consistent |
| `radomski1989_extraction.md` (sugar-ring C-2 relaxation in minutes, pH 3.7-4.0, 37 °C) | ~0.2 min⁻¹ (glyceraldehyde ring) | 2.4 × 10⁻² min⁻¹ on the acid limb, 30 °C | same order on the acid side; the penicilloate's electron-withdrawing C-2 substituent slows it |
| `kallen1971b_extraction.md` (TC k_r rises with pH from 2.8 to 5.6, flat to 6.6, falls above) | formaldehyde ring release | aliphatic ring opening falls from pH 3 to 6, flat to 12 | different steps (release vs opening); not comparable directly |

## 7. Flags

1. **[!] Epimerisation is a lower bound on ring opening** (if closure is faster than conformational change, the
   authors note); the reported numbers are epimerisation constants.
2. **[!] Buffers unnamed, concentrations unprinted, no errors**; no buffer-dependence test.
3. 30 °C and I 1.0 M only.
4. Penicillin numbering: C-5 here = thiazolidine C-2; C-3 = thiazolidine C-4 (carboxylate); C-2 = thiazolidine C-5
   (gem-dimethyl).

## M1c role (frozen bd641ae8)

**Role fixed by the case (sec. 4.3):** band source for the pH shape of ring opening (R-rel).

**Band it can supply:**
- **Coordinate: R-rel's n_r** (the release's pH slope; case band 0-1, analogues arguing ≤ 0). For an aliphatic
  thiazolidine-4-carboxylate at 30 °C, I 1.0 M: ring OPENING has slope **0 over pH 6-12** and **−0.83 dex across
  the ring-N pKa (≈ 4.6-5.1)**, i.e. negative between pH ~4 and ~6; the hydroxide term (n = 1) enters only above pH
  ~12.7. Hydrolytic RELEASE is ≥ 20x slower than opening below pH 7 and undetected above pH 8: negative slope.
  Both say n_r ≤ 0 over the pots' pH 4-8, which is the analogues' verdict (Karpinski, Bustillo Trueba, Nagasawa,
  Pesek) from a fifth, independent laboratory.
- **Secondary, P1a / P2 exchange speed:** opening ≥ 3.0 × 10⁻³ min⁻¹ (neutral ring) and 2.4 × 10⁻² min⁻¹ (protonated)
  at 30 °C for a slow (gem-dimethyl, electron-poor C-2) member of the class, consistent with the case's ≥ 2 × 10⁻³
  min⁻¹ band (Nagasawa) and faster at process temperature.
- **Transfer:** same reaction class (thiazolidine-4-carboxylate, aliphatic C-2), different substituents; 30 °C. The
  pH shape transfers better than any magnitude.
- **Not decisive.** 30 °C, no cysteine-pentose ring, no release rate. Nothing the case would call decisive.
